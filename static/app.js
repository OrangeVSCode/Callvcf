const state = { metadata: null, activeTab: "existence", lastResult: null };
const $ = (id) => document.getElementById(id);
const categoryOrder = ["HOM_REF", "HET", "HOM_ALT", "MISSING", "OTHER"];
const shortLabels = { HOM_REF: "0/0", HET: "0/1", HOM_ALT: "1/1", MISSING: "缺失", OTHER: "其他" };

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>'"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","'":"&#39;",'"':"&quot;"}[c]));
}

function formatBytes(n) {
  if (n == null) return "—";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let value = Number(n), i = 0;
  while (value >= 1024 && i < units.length - 1) { value /= 1024; i++; }
  return `${value >= 10 || i === 0 ? value.toFixed(0) : value.toFixed(1)} ${units[i]}`;
}

function formatNumber(n) { return n == null ? "—" : Number(n).toLocaleString("zh-CN"); }

function toast(message, error = false) {
  const el = $("toast");
  el.textContent = message;
  el.className = `toast show${error ? " error" : ""}`;
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => { el.className = "toast"; }, 3800);
}

async function api(route, payload) {
  const response = await fetch(route, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const data = await response.json().catch(() => ({ ok: false, error: `HTTP ${response.status}` }));
  if (!response.ok || !data.ok) throw new Error(data.error || `HTTP ${response.status}`);
  return data.data;
}

function currentPath() {
  const path = $("vcfPath").value.trim();
  if (!path) throw new Error("请先输入并载入 VCF 路径");
  return path;
}

function lociText(required = true) {
  const text = $("lociInput").value.trim();
  if (required && !text) throw new Error("请输入至少一个位点");
  return text;
}

function selectedSamples(required = true) {
  const values = $("samplesInput").value.split(/[\s,;]+/).map(x => x.trim()).filter(Boolean);
  const unique = [...new Set(values)];
  if (required && !unique.length) throw new Error("请输入至少一个样本");
  return unique;
}

function setLoading(target, button, message = "正在查询…") {
  target.className = "result-body loading";
  target.textContent = message;
  button.disabled = true;
  const old = button.textContent;
  button.dataset.oldText = old;
  button.textContent = "处理中…";
}

function clearLoading(target, button) {
  target.className = "result-body";
  button.disabled = false;
  button.textContent = button.dataset.oldText || button.textContent;
}

async function checkHealth() {
  const status = $("serverStatus");
  try {
    const response = await fetch("/api/health");
    if (!response.ok) throw new Error();
    const data = await response.json();
    status.className = "status-pill online";
    status.innerHTML = `<span></span>本地服务 · ${escapeHtml(data.backend || "ready")}`;
  } catch {
    status.className = "status-pill error";
    status.innerHTML = "<span></span>连接失败";
  }
}

function renderMetadata(meta) {
  const types = Object.entries(meta.observed_types || {}).map(([k, v]) => `<span class="chip">${escapeHtml(k)} · ${formatNumber(v)}</span>`).join("");
  $("metadata").className = "metadata";
  $("metadata").innerHTML = `
    <div class="meta-grid">
      <div class="meta-card"><span>存储格式</span><strong title="${escapeHtml(meta.storage)}">${escapeHtml(meta.storage)}</strong></div>
      <div class="meta-card"><span>VCF 版本</span><strong>${escapeHtml(meta.vcf_version)}</strong></div>
      <div class="meta-card"><span>样本数</span><strong>${formatNumber(meta.sample_count)}</strong></div>
      <div class="meta-card"><span>记录数</span><strong>${formatNumber(meta.record_count)}</strong></div>
      <div class="meta-card"><span>染色体/Contig</span><strong>${formatNumber(meta.contig_count)}</strong></div>
      <div class="meta-card"><span>查询模式</span><strong title="${escapeHtml(meta.query_mode)}">${meta.indexed ? "已索引 · 随机访问" : "未索引 · 流式筛选"}</strong></div>
    </div>
    <div class="meta-path">${escapeHtml(meta.path)} · ${formatBytes(meta.file_size)}</div>
    <div class="type-chips">${types || '<span class="chip">未观察到记录</span>'}</div>`;
  $("samplePicker").classList.toggle("hidden", !meta.samples.length);
  renderSampleSuggestions("");
}

async function inspectVCF() {
  const button = $("inspectBtn");
  button.disabled = true; button.textContent = "识别中…";
  $("metadata").className = "metadata loading"; $("metadata").textContent = "";
  try {
    const meta = await api("/api/inspect", { path: currentPath() });
    state.metadata = meta;
    localStorage.setItem("vcfExplorerPath", meta.path);
    renderMetadata(meta);
    toast(`已载入 ${meta.name}，${meta.sample_count} 个样本`);
  } catch (error) {
    $("metadata").className = "metadata empty-state";
    $("metadata").textContent = error.message;
    toast(error.message, true);
  } finally {
    button.disabled = false; button.textContent = "载入并识别";
  }
}

async function selectLocalFile() {
  const button = $("selectFileBtn");
  button.disabled = true; button.textContent = "等待选择…";
  try {
    const result = await api("/api/select-file", {});
    if (result.path) {
      $("vcfPath").value = result.path;
      await inspectVCF();
    }
  } catch (error) { toast(error.message, true); }
  finally { button.disabled = false; button.textContent = "选择文件…"; }
}

async function shutdownLocal() {
  if (!confirm("关闭本地 CallVCF 工具？")) return;
  try { await api("/api/shutdown", {}); } catch (_) {}
  document.body.innerHTML = '<main style="max-width:720px;margin:12vh auto;padding:30px"><section class="panel"><h1 style="font-size:42px">CallVCF 已关闭</h1><p class="lede">可以安全关闭这个页面。下次双击 CallVCF 即可重新启动。</p></section></main>';
}

async function discoverFiles() {
  const button = $("discoverBtn");
  button.disabled = true; button.textContent = "扫描中…";
  try {
    const result = await api("/api/files", { root: $("discoverRoot").value.trim(), max_depth: 5 });
    const select = $("fileSelect");
    select.innerHTML = result.files.map(f => `<option value="${escapeHtml(f.path)}">${escapeHtml(f.path)} · ${formatBytes(f.size)}</option>`).join("");
    $("fileSelectWrap").classList.toggle("hidden", !result.files.length);
    if (result.files.length) {
      $("vcfPath").value = result.files[0].path;
      toast(`发现 ${result.files.length} 个 VCF/BCF 文件${result.truncated ? "（结果已截断）" : ""}`);
    } else toast("目录中未发现 VCF/BCF 文件", true);
  } catch (error) { toast(error.message, true); }
  finally { button.disabled = false; button.textContent = "扫描 VCF 文件"; }
}

function renderSampleSuggestions(query) {
  if (!state.metadata) return;
  const selected = new Set(selectedSamples(false));
  const q = query.trim().toLowerCase();
  const matches = state.metadata.samples.filter(s => !q || s.toLowerCase().includes(q)).slice(0, 200);
  $("sampleSearchCount").textContent = `${matches.length}${state.metadata.sample_count > 200 ? " / 最多显示 200" : ""}`;
  $("sampleSuggestions").innerHTML = matches.map(s => `<button class="sample-token${selected.has(s) ? " added" : ""}" data-sample="${escapeHtml(s)}">${escapeHtml(s)}</button>`).join("");
}

function toggleSample(sample) {
  const samples = selectedSamples(false);
  const index = samples.indexOf(sample);
  if (index >= 0) samples.splice(index, 1); else samples.push(sample);
  $("samplesInput").value = samples.join("\n");
  renderSampleSuggestions($("sampleSearch").value);
}

function renderExistence(data) {
  const rows = data.results.flatMap(item => item.exists
    ? item.records.map((r, i) => `<tr><td>${escapeHtml(item.query)}</td><td class="yes">存在${item.records.length > 1 ? ` · 记录 ${i + 1}` : ""}</td><td>${escapeHtml(r.ref)}</td><td>${escapeHtml(r.alt)}</td><td><span class="variant-badge">${escapeHtml(r.variant_type)}${r.svtype ? ` · ${escapeHtml(r.svtype)}` : ""}</span></td><td>${escapeHtml(r.id || ".")}</td><td>${escapeHtml(r.end)}</td></tr>`)
    : [`<tr><td>${escapeHtml(item.query)}</td><td class="no">不存在</td><td colspan="5">—</td></tr>`]
  ).join("");
  return `<div class="table-wrap"><table><thead><tr><th>查询位点</th><th>结果</th><th>REF</th><th>ALT</th><th>类型</th><th>ID</th><th>END</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}

function distributionCard(item) {
  const total = item.total_samples || 1;
  const segments = categoryOrder.map(k => `<span class="${k}" style="width:${(item.counts[k] || 0) / total * 100}%" title="${shortLabels[k]}: ${item.counts[k] || 0}"></span>`).join("");
  const counts = categoryOrder.map(k => `<div class="count-item"><strong>${formatNumber(item.counts[k] || 0)}</strong><span>${shortLabels[k]}</span></div>`).join("");
  const sampleRows = categoryOrder.map(k => `<div class="sample-list-row"><b>${shortLabels[k]} (${item.counts[k] || 0})</b><span>${(item.samples[k] || []).map(escapeHtml).join(", ") || "—"}</span></div>`).join("");
  const r = item.record;
  return `<article class="result-card">
    <div class="result-card-head"><div><h4>${escapeHtml(r.key)}</h4><p class="sub">${escapeHtml(r.variant_type)}${r.svtype ? ` / ${escapeHtml(r.svtype)}` : ""} · ${escapeHtml(r.id || "无 ID")}</p></div><span class="variant-badge">n=${item.total_samples}</span></div>
    <div class="result-card-body"><div class="stacked-bar">${segments}</div><div class="count-grid">${counts}</div></div>
    <details><summary>查看各类别的具体样本</summary><div class="sample-lists">${sampleRows}</div></details>
  </article>`;
}

function renderDistribution(data) {
  const missing = data.missing_loci.length ? `<div class="missing-box">VCF 中未找到：${data.missing_loci.map(escapeHtml).join("、")}</div>` : "";
  const cards = data.records.map(distributionCard).join("");
  return `${missing}${cards ? `<div class="distribution-grid">${cards}</div>` : '<div class="empty-state">没有匹配的变异记录。</div>'}`;
}

function renderMatrix(data) {
  const missing = data.missing_loci.length ? `<div class="missing-box">VCF 中未找到：${data.missing_loci.map(escapeHtml).join("、")}</div>` : "";
  if (!data.columns.length) return `${missing}<div class="empty-state">没有可显示的变异记录。</div>`;
  const head = data.columns.map(c => `<th title="${escapeHtml(c.ref)}>${escapeHtml(c.alt)} · ${escapeHtml(c.variant_type)}">${escapeHtml(c.key)}</th>`).join("");
  const rows = data.rows.map(row => `<tr><th>${escapeHtml(row.sample)}</th>${row.values.map(v => `<td class="gt-cell gt-${v.category}" title="${escapeHtml(v.label)}">${escapeHtml(v.gt)} · ${escapeHtml(v.alleles)}</td>`).join("")}</tr>`).join("");
  return `${missing}<div class="table-wrap"><table><thead><tr><th>样本</th>${head}</tr></thead><tbody>${rows}</tbody></table></div>`;
}

function statsCard(item) {
  const total = item.total_records || 1;
  const segments = categoryOrder.map(k => `<span class="${k}" style="width:${(item.counts[k] || 0) / total * 100}%" title="${shortLabels[k]}: ${item.counts[k] || 0}"></span>`).join("");
  const counts = categoryOrder.map(k => `<div class="count-item"><strong>${formatNumber(item.counts[k] || 0)}</strong><span>${shortLabels[k]}</span></div>`).join("");
  const typeRows = Object.entries(item.by_variant_type || {}).map(([type, c]) => `<tr><th>${escapeHtml(type)}</th>${categoryOrder.map(k => `<td>${formatNumber(c[k] || 0)}</td>`).join("")}</tr>`).join("");
  return `<article class="result-card">
    <div class="result-card-head"><div><h4>${escapeHtml(item.sample)}</h4><p class="sub">共统计 ${formatNumber(item.total_records)} 条记录</p></div></div>
    <div class="result-card-body"><div class="stacked-bar">${segments}</div><div class="count-grid">${counts}</div>
      <div class="stat-types"><table><thead><tr><th>变异类型</th>${categoryOrder.map(k => `<th>${shortLabels[k]}</th>`).join("")}</tr></thead><tbody>${typeRows}</tbody></table></div>
    </div>
  </article>`;
}

function genericTable(rows, preferred = []) {
  if (!rows || !rows.length) return '<div class="empty-state">没有匹配记录。</div>';
  const keys = [...new Set(rows.flatMap(row => Object.keys(row)))];
  const ordered = [...preferred.filter(k => keys.includes(k)), ...keys.filter(k => !preferred.includes(k))];
  const head = ordered.map(k => `<th>${escapeHtml(k)}</th>`).join("");
  const body = rows.map(row => `<tr>${ordered.map(k => {
    const value = typeof row[k] === "object" && row[k] !== null ? JSON.stringify(row[k]) : row[k];
    return `<td>${escapeHtml(value ?? "")}</td>`;
  }).join("")}</tr>`).join("");
  return `<div class="table-wrap"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`;
}

function renderLdChart(ld) {
  const source = ld.points || [];
  if (!source.length) return '<div class="empty-state">没有可绘制的 LD 点。</div>';
  const maxPoints = 2200;
  let points = source;
  if (source.length > maxPoints) {
    const linkedKeys = new Set((ld.linked_variants || []).map(x => x.key));
    const step = Math.ceil(source.length / maxPoints);
    points = source.filter((x, i) => i % step === 0 || linkedKeys.has(x.key));
  }
  const region = ld.search_region;
  const left = 54, right = 18, top = 12, bottom = 34, width = 900, height = 210;
  const plotW = width - left - right, plotH = height - top - bottom;
  const span = Math.max(1, region.end - region.start);
  const x = pos => left + (pos - region.start) / span * plotW;
  const y = r2 => top + (1 - r2) * plotH;
  const linkedKeys = new Set((ld.linked_variants || []).map(item => item.key));
  const circles = points.map(item => {
    const cls = item.pos === ld.lead_record.pos ? "lead-point" : linkedKeys.has(item.key) ? "linked" : "point";
    const radius = item.pos === ld.lead_record.pos ? 5 : 2.6;
    return `<circle class="${cls}" cx="${x(item.pos).toFixed(2)}" cy="${y(item.r2).toFixed(2)}" r="${radius}"><title>${escapeHtml(item.key)} · r²=${item.r2} · n=${item.n}</title></circle>`;
  }).join("");
  const ticks = [0, .25, .5, .75, 1].map(v => `<line class="grid" x1="${left}" x2="${width-right}" y1="${y(v)}" y2="${y(v)}"></line><text x="${left-9}" y="${y(v)+3}" text-anchor="end">${v}</text>`).join("");
  const thresholdY = y(ld.threshold);
  return `<div class="ld-chart"><svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Lead 位点与窗口内变异的 r² 散点图">
    ${ticks}<line class="axis" x1="${left}" x2="${width-right}" y1="${height-bottom}" y2="${height-bottom}"></line>
    <line x1="${left}" x2="${width-right}" y1="${thresholdY}" y2="${thresholdY}" stroke="#c9494d" stroke-dasharray="5 4"></line>
    <text x="${width-right}" y="${thresholdY-5}" text-anchor="end">阈值 r²=${ld.threshold}</text>
    ${circles}<text x="${left}" y="${height-9}">${escapeHtml(region.chrom)}:${formatNumber(region.start)}</text>
    <text x="${width-right}" y="${height-9}" text-anchor="end">${escapeHtml(region.chrom)}:${formatNumber(region.end)}</text>
    <text x="${width/2}" y="${height-9}" text-anchor="middle">基因组位置</text>
  </svg></div>`;
}

function renderLeadResult(data) {
  const lead = data.lead_record || {};
  const ld = data.ld;
  const regionText = ld ? `${ld.linkage_region.chrom}:${formatNumber(ld.linkage_region.start)}–${formatNumber(ld.linkage_region.end)}` : "未计算";
  let html = `<div class="lead-summary">
    <div class="meta-card"><span>Lead 记录</span><strong>${escapeHtml(lead.key || data.lead_locus)}</strong></div>
    <div class="meta-card"><span>变异类型</span><strong>${escapeHtml(lead.variant_type || "—")}${lead.svtype ? ` / ${escapeHtml(lead.svtype)}` : ""}</strong></div>
    <div class="meta-card"><span>连锁区域</span><strong>${regionText}</strong></div>
    <div class="meta-card"><span>连锁变异 / SNP</span><strong>${ld ? `${formatNumber(ld.linked_variant_count)} / ${formatNumber(ld.linked_snp_count)}` : "未计算"}</strong></div>
  </div>`;

  if (ld) {
    const linkedRows = (ld.linked_variants || []).map(x => ({位点: `${x.chrom}:${x.pos}`, ID: x.id || ".", REF: x.ref, ALT: x.alt, 类型: x.variant_type, "r²": x.r2, 有效样本: x.n, 距Lead_bp: x.distance}));
    html += `<section class="analysis-section"><h4>LD 区段与连锁变异集</h4><p class="sub">搜索 ${ld.search_region.chrom}:${formatNumber(ld.search_region.start)}–${formatNumber(ld.search_region.end)}；检测 ${formatNumber(ld.tested_variant_count)} 个二等位变异。绿色为达到阈值的变异，红色为 Lead。</p>
      ${renderLdChart(ld)}
      <div class="pane-intro" style="margin:12px 0 8px"><p>连锁跨度 ${formatNumber(ld.linkage_region.span_bp)} bp；r² ≥ ${ld.threshold}。</p><div><button class="secondary" id="copyLinkedBtn">复制位点</button> <button class="secondary" id="downloadLinkedBtn">下载 TSV</button></div></div>
      ${genericTable(linkedRows, ["位点", "ID", "REF", "ALT", "类型", "r²", "有效样本", "距Lead_bp"])}</section>`;
  }

  if (data.gene) {
    const rows = (data.gene.matches || []).map(x => ({基因: x.name, 基因ID: x.id, 关系: x.relation, 距离_bp: x.distance, 染色体起点: x.start, 染色体终点: x.end, 链: x.strand, 重叠功能元件: x.overlapping_feature?.type || ""}));
    html += `<section class="analysis-section"><h4>相对基因位置</h4><p class="sub">${data.gene.status === "not_configured" ? "未配置 GFF3/GTF。" : escapeHtml(data.gene.path || "")}</p>${genericTable(rows, ["基因", "基因ID", "关系", "距离_bp", "染色体起点", "染色体终点", "链", "重叠功能元件"])}</section>`;
  }

  if (data.function) {
    html += `<section class="analysis-section"><h4>突变功能</h4><div class="result-split"><div><p class="sub">VCF INFO（ANN / CSQ / BCSQ）</p>${genericTable(data.function.vcf_info, ["source", "effect", "impact", "gene", "gene_id", "feature", "hgvs_c", "hgvs_p"])}</div><div><p class="sub">外部功能注释表</p>${genericTable(data.function.external_table)}</div></div></section>`;
  }
  if (data.domain) {
    html += `<section class="analysis-section"><h4>蛋白结构域</h4><p class="sub">用于匹配的基因/蛋白标识：${(data.domain.identifiers || []).map(escapeHtml).join("、") || "未从注释中获得标识"}</p>${genericTable(data.domain.matches)}</section>`;
  }
  if (data.phenotype) {
    const traits = (data.phenotype.traits || []).map(x => ({表型: x.trait, 区域记录数: x.record_count, 最小P值: x.min_p, 峰值Marker: x.top_marker, 峰值位置: x.top_pos, 文件: x.source_file}));
    html += `<section class="analysis-section"><h4>区域表型关联</h4><p class="sub">读取 ${formatNumber((data.phenotype.files || []).length)} 个 .ps 文件，区域内匹配 ${formatNumber((data.phenotype.records || []).length)} 条记录。</p>${genericTable(traits, ["表型", "区域记录数", "最小P值", "峰值Marker", "峰值位置", "文件"])}</section>`;
  }
  if (data.ldblockshow) {
    const run = data.ldblockshow;
    html += `<section class="analysis-section"><h4>LDBlockShow</h4><p class="${run.ok ? "status-good" : "status-bad"}">${run.ok ? "运行完成" : `未完成：${escapeHtml(run.error || `返回码 ${run.return_code}`)}`}</p>
      ${run.command ? `<pre class="code-block">${escapeHtml(run.command.join(" "))}</pre>` : ""}
      ${run.output_files?.length ? genericTable(run.output_files.map(x => ({输出文件: x}))) : ""}
      ${run.stderr ? `<details><summary>查看软件日志</summary><pre class="code-block">${escapeHtml(run.stderr)}</pre></details>` : ""}</section>`;
  }
  return html;
}

function saveAdvancedSettings() {
  const ids = ["gffPath", "annotationPath", "domainPath", "phenotypePath", "ldblockshowPath", "outputDir", "ldWindow", "ldThreshold", "ldMinSamples"];
  localStorage.setItem("vcfExplorerAdvanced", JSON.stringify(Object.fromEntries(ids.map(id => [id, $(id).value]))));
}

function restoreAdvancedSettings() {
  try {
    const saved = JSON.parse(localStorage.getItem("vcfExplorerAdvanced") || "{}");
    Object.entries(saved).forEach(([id, value]) => { if ($(id)) $(id).value = value; });
  } catch (_) {}
}

async function browseResource(button) {
  button.disabled = true;
  try {
    const result = await api("/api/select-resource", { kind: button.dataset.kind });
    if (result.path) { $(button.dataset.target).value = result.path; saveAdvancedSettings(); }
  } catch (error) { toast(error.message, true); }
  finally { button.disabled = false; }
}

async function runLeadAnalysis() {
  const button = $("leadRunBtn"), target = $("leadResult");
  try {
    const leadLocus = $("leadLocus").value.trim() || lociText(false).split(/[\s,;]+/).find(Boolean);
    if (!leadLocus) throw new Error("请输入一个 Lead 位点");
    const options = { ld: $("optLd").checked, gene: $("optGene").checked, function: $("optFunction").checked,
      domain: $("optDomain").checked, phenotype: $("optPhenotype").checked, ldblockshow: $("optLdblockshow").checked };
    if (!Object.values(options).some(Boolean)) throw new Error("请至少勾选一项分析功能");
    const payload = { path: currentPath(), lead_locus: leadLocus, window_kb: $("ldWindow").value,
      r2_threshold: $("ldThreshold").value, min_samples: $("ldMinSamples").value, options,
      gff_path: $("gffPath").value.trim(), annotation_path: $("annotationPath").value.trim(),
      domain_path: $("domainPath").value.trim(), phenotype_path: $("phenotypePath").value.trim(),
      ldblockshow_path: $("ldblockshowPath").value.trim(), ldblockshow_wsl: $("ldblockshowWsl").checked,
      output_dir: $("outputDir").value.trim() };
    saveAdvancedSettings();
    setLoading(target, button, options.ldblockshow ? "正在计算 LD 并运行 LDBlockShow…" : "正在分析 Lead 位点…");
    const result = await api("/api/lead-analysis", payload);
    state.lastLeadResult = result;
    clearLoading(target, button);
    target.innerHTML = renderLeadResult(result);
    toast("Lead 位点分析完成");
  } catch (error) {
    clearLoading(target, button);
    target.className = "result-body empty-state";
    target.textContent = error.message;
    toast(error.message, true);
  }
}

function linkedTsv(data) {
  const rows = data?.ld?.linked_variants || [];
  return ["chrom\tpos\tid\tref\talt\ttype\tr2\tn\tdistance_bp", ...rows.map(x => [x.chrom, x.pos, x.id || ".", x.ref, x.alt, x.variant_type, x.r2, x.n, x.distance].join("\t"))].join("\n");
}

async function copyLinked() {
  const text = (state.lastLeadResult?.ld?.linked_variants || []).map(x => `${x.chrom}:${x.pos}`).join("\n");
  await navigator.clipboard.writeText(text);
  toast(`已复制 ${text ? text.split("\n").length : 0} 个位点`);
}

function downloadLinked() {
  const blob = new Blob([linkedTsv(state.lastLeadResult)], {type: "text/tab-separated-values;charset=utf-8"});
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `linked_variants_${(state.lastLeadResult?.lead_locus || "lead").replace(":", "_")}.tsv`;
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 1000);
}

async function runAction(action, button) {
  const map = {
    check: ["existenceResult", "/api/check"],
    distribution: ["distributionResult", "/api/distribution"],
    matrix: ["matrixResult", "/api/sample-loci"],
    stats: ["statsResult", "/api/sample-stats"],
  };
  const [targetId, route] = map[action];
  const target = $(targetId);
  try {
    const payload = { path: currentPath() };
    if (action !== "stats" || lociText(false)) payload.loci = lociText(action !== "stats");
    if (["matrix", "stats"].includes(action)) payload.samples = selectedSamples(true);
    setLoading(target, button, action === "stats" && !payload.loci ? "正在扫描整个 VCF…" : "正在查询…");
    const result = await api(route, payload);
    state.lastResult = result;
    clearLoading(target, button);
    if (action === "check") target.innerHTML = renderExistence(result);
    if (action === "distribution") target.innerHTML = renderDistribution(result);
    if (action === "matrix") target.innerHTML = renderMatrix(result);
    if (action === "stats") target.innerHTML = `<p class="sub">处理记录：${formatNumber(result.processed_records)}</p><div class="stats-grid">${result.results.map(statsCard).join("")}</div>`;
    toast("查询完成");
  } catch (error) {
    clearLoading(target, button);
    target.className = "result-body empty-state";
    target.textContent = error.message;
    toast(error.message, true);
  }
}

function initEvents() {
  $("selectFileBtn").addEventListener("click", selectLocalFile);
  $("shutdownBtn").addEventListener("click", shutdownLocal);
  $("inspectBtn").addEventListener("click", inspectVCF);
  $("discoverBtn").addEventListener("click", discoverFiles);
  $("fileSelect").addEventListener("change", e => { $("vcfPath").value = e.target.value; });
  $("sampleSearch").addEventListener("input", e => renderSampleSuggestions(e.target.value));
  $("samplesInput").addEventListener("input", () => renderSampleSuggestions($("sampleSearch").value));
  $("sampleSuggestions").addEventListener("click", e => { const btn = e.target.closest("[data-sample]"); if (btn) toggleSample(btn.dataset.sample); });
  $("leadRunBtn").addEventListener("click", runLeadAnalysis);
  document.querySelectorAll(".browse-resource").forEach(btn => btn.addEventListener("click", () => browseResource(btn)));
  $("leadResult").addEventListener("click", e => {
    if (e.target.closest("#copyLinkedBtn")) copyLinked().catch(err => toast(err.message, true));
    if (e.target.closest("#downloadLinkedBtn")) downloadLinked();
  });
  document.querySelectorAll(".tab").forEach(tab => tab.addEventListener("click", () => {
    document.querySelectorAll(".tab").forEach(x => x.classList.toggle("active", x === tab));
    document.querySelectorAll(".tab-pane").forEach(x => x.classList.toggle("active", x.id === `pane-${tab.dataset.tab}`));
    state.activeTab = tab.dataset.tab;
  }));
  document.querySelectorAll(".run-btn").forEach(btn => btn.addEventListener("click", () => runAction(btn.dataset.action, btn)));
}

document.addEventListener("DOMContentLoaded", () => {
  const saved = localStorage.getItem("vcfExplorerPath");
  if (saved) $("vcfPath").value = saved;
  initEvents();
  restoreAdvancedSettings();
  checkHealth();
});
