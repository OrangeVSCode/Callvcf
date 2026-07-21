const state = { metadata: null, activeTab: "existence", lastResult: null, lastLeadResult: null, lastProfileResult: null, qualityCatalog: null, qualityRun: null, qualityPoll: null, repairCatalog: null, repairPlan: null, repairPoll: null };
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
function formatPercent(n) { return n == null ? "—" : `${(Number(n) * 100).toFixed(2)}%`; }

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
  const queryLabel = meta.index_usable
    ? "已索引 · 随机访问"
    : (meta.indexed ? "检测到索引 · 当前流式读取" : (meta.compressed ? "压缩直读 · 流式筛选" : "未索引 · 流式筛选"));
  const recordCount = meta.record_count == null
    ? `<strong>—</strong><button class="meta-action" id="countRecordsBtn" title="顺序读取压缩数据流，不生成解压文件">统计总数</button>`
    : `<strong>${formatNumber(meta.record_count)}</strong><small>${meta.record_count_source === "index" ? "来自索引" : "完整流式统计"}</small>`;
  $("metadata").className = "metadata";
  $("metadata").innerHTML = `
    <div class="meta-grid">
      <div class="meta-card"><span>存储格式</span><strong title="${escapeHtml(meta.storage)}">${escapeHtml(meta.storage)}</strong></div>
      <div class="meta-card"><span>VCF 版本</span><strong>${escapeHtml(meta.vcf_version)}</strong></div>
      <div class="meta-card"><span>样本数</span><strong>${formatNumber(meta.sample_count)}</strong></div>
      <div class="meta-card record-count-card"><span>记录数</span>${recordCount}</div>
      <div class="meta-card"><span>染色体/Contig</span><strong>${formatNumber(meta.contig_count)}</strong></div>
      <div class="meta-card"><span>查询模式</span><strong title="${escapeHtml(meta.query_mode)}">${queryLabel}</strong></div>
      <div class="meta-card"><span>磁盘占用</span><strong title="${escapeHtml(meta.space_mode || "")}">原文件直读 · 不生成解压副本</strong></div>
    </div>
    <div class="meta-path">${escapeHtml(meta.path)} · ${formatBytes(meta.file_size)}</div>
    <div class="type-chips">${types || '<span class="chip">未观察到记录</span>'}</div>`;
  $("samplePicker").classList.toggle("hidden", !meta.samples.length);
  renderSampleSuggestions("");
  const countButton = $("countRecordsBtn");
  if (countButton) countButton.addEventListener("click", countRecords);
}

async function countRecords() {
  const button = $("countRecordsBtn");
  if (!button || !state.metadata) return;
  button.disabled = true;
  button.textContent = "统计中…";
  toast("正在完整读取压缩数据流统计记录数；不会生成解压文件");
  try {
    const result = await api("/api/count-records", { path: currentPath() });
    state.metadata.record_count = result.record_count;
    state.metadata.record_count_source = result.method === "index" ? "index" : "full_stream";
    renderMetadata(state.metadata);
    const elapsed = result.cached ? "已使用缓存" : `耗时 ${result.elapsed_seconds} 秒`;
    toast(`记录总数：${formatNumber(result.record_count)}（${elapsed}）`);
  } catch (error) {
    button.disabled = false;
    button.textContent = "重新统计";
    toast(error.message, true);
  }
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

function r2Color(value) {
  if (value == null) return "rgb(238,241,239)";
  const v = Math.max(0, Math.min(1, Number(value)));
  const r = Math.round(255 - 18 * (1 - v));
  const g = Math.round(246 - 202 * v);
  const b = Math.round(220 - 184 * v);
  return `rgb(${r},${g},${b})`;
}

function geneTrackMarkup(data, x, top, plotW) {
  const models = (data.gene_track?.models || []).slice(0, 24);
  if (!models.length) return {height: 44, markup: `<text x="${x}" y="${top+22}" fill="#71827a">Gene structure: provide a GFF3/GTF file and select “基因结构轨道”</text>`};
  const region = data.ld.search_region, span = Math.max(1, region.end-region.start);
  const px = pos => x+(pos-region.start)/span*plotW;
  const lanes = [];
  models.forEach(model => {
    let lane = lanes.findIndex(end => model.start > end);
    if (lane < 0) { lane = lanes.length; lanes.push(model.end); } else lanes[lane] = model.end;
    model._lane = Math.min(lane, 5);
  });
  const height = 50 + Math.min(6, lanes.length)*30;
  const markup = models.filter(model => model._lane < 6).map(model => {
    const y = top+42+model._lane*30, x1=px(model.start), x2=px(model.end);
    const direction = model.strand === "-" ? -1 : 1;
    const arrows=[];
    for(let ax=x1+18; ax<x2-8; ax+=32) arrows.push(`<path d="M ${ax} ${y-3} l ${6*direction} 3 l ${-6*direction} 3" fill="none" stroke="#46635a"/>`);
    const features=(model.features||[]).filter(f=>["exon","cds","five_prime_utr","three_prime_utr","utr"].includes(f.type)).map(f=>{
      const fx1=px(f.start), fx2=px(f.end), isCds=f.type==="cds", h=isCds?14:9;
      return `<rect x="${Math.min(fx1,fx2)}" y="${y-h/2}" width="${Math.max(1,Math.abs(fx2-fx1))}" height="${h}" fill="${isCds?'#2c6b5a':'#8fb9aa'}"><title>${escapeHtml(model.name)} · ${f.type} · ${f.start}-${f.end}</title></rect>`;
    }).join("");
    return `<g><line x1="${x1}" x2="${x2}" y1="${y}" y2="${y}" stroke="#46635a" stroke-width="1.5"/>${arrows.join("")}${features}<text x="${x1}" y="${y-10}" fill="#344b43">${escapeHtml(model.name)} (${escapeHtml(model.strand)})</text></g>`;
  }).join("");
  return {height, markup:`<text x="${x}" y="${top+12}" fill="#344b43">Gene structure</text>${markup}`};
}

function renderRegionalLdMetric(data, metric) {
  const ld = data.ld, heatmap = ld?.heatmap;
  if (!heatmap?.variants?.length) return '<div class="empty-state">达到阈值的位点不足，无法绘制 LD 热图。</div>';
  const width = 1000, left = 76, right = 24, plotW = width-left-right;
  const region = ld.search_region, span = Math.max(1, region.end - region.start);
  const x = pos => left + (pos - region.start) / span * plotW;
  const isDprime = metric === "dprime", matrix = isDprime ? heatmap.matrix_dprime : (heatmap.matrix_r2 || heatmap.matrix);
  const metricLabel = isDprime ? "D′" : "R²", metricSlug = isDprime ? "Dprime" : "R2";
  const association = (data.phenotype?.records || []).filter(r => r.chrom === region.chrom && r.pos >= region.start && r.pos <= region.end);
  const maxLogP = Math.max(1, ...association.map(r => -Math.log10(Math.max(Number(r.pvalue), 1e-300))));
  const assocTop = 24, assocBottom = 185, assocH = assocBottom - assocTop;
  const assocMarks = association.slice(0, 12000).map(r => {
    const value = -Math.log10(Math.max(Number(r.pvalue), 1e-300));
    return `<circle cx="${x(r.pos).toFixed(2)}" cy="${(assocBottom - value / maxLogP * assocH).toFixed(2)}" r="2.2" fill="#7d9088"><title>${escapeHtml(r.trait)} · ${escapeHtml(r.marker)} · P=${r.pvalue}</title></circle>`;
  }).join("");
  const geneTrack = geneTrackMarkup(data, left, 218, plotW), blockTop=218+geneTrack.height;
  const leads = ld.lead_results?.length ? ld.lead_results : [{lead_record: ld.lead_record, block: ld.linkage_region}];
  const blockMarks = leads.map((item, i) => {
    const y = blockTop+18+(i%3)*13, block=item.block||ld.linkage_region;
    return `<line x1="${x(block.start)}" x2="${x(block.end)}" y1="${y}" y2="${y}" stroke="#d14c4c" stroke-width="5" opacity=".72"><title>${escapeHtml(item.lead_record.key)} · ${block.start}-${block.end}</title></line><path d="M ${x(item.lead_record.pos)-5} ${y-9} L ${x(item.lead_record.pos)+5} ${y-9} L ${x(item.lead_record.pos)} ${y-1} Z" fill="#9f252d"/>`;
  }).join("");
  const variants=heatmap.variants, n=variants.length, heatTop=blockTop+76, heatH=300;
  const height=heatTop+heatH+76, unitX=plotW/Math.max(1,n), unitY=heatH/Math.max(1,n);
  const cells = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      const value = matrix[i]?.[j];
      const cx = left + (i + j + 1) * unitX / 2;
      const cy = heatTop + (i - j + 1) * unitY / 2;
      const dx = unitX / 2 + .35, dy = unitY / 2 + .35;
      cells.push(`<path d="M ${cx} ${cy-dy} L ${cx+dx} ${cy} L ${cx} ${cy+dy} L ${cx-dx} ${cy} Z" fill="${r2Color(value)}"><title>${escapeHtml(variants[j].key)} × ${escapeHtml(variants[i].key)} · ${metricLabel}=${value == null ? "NA" : value}</title></path>`);
    }
  }
  const leadGuides=leads.map(item=>`<line x1="${x(item.lead_record.pos)}" x2="${x(item.lead_record.pos)}" y1="${assocTop}" y2="${heatTop-18}" stroke="#c64048" stroke-width="1.2" stroke-dasharray="6 5"/>`).join("");
  const yTicks = [0, .25, .5, .75, 1].map(v => `<line x1="${left}" x2="${width-right}" y1="${assocBottom-v*assocH}" y2="${assocBottom-v*assocH}" stroke="#e0e7e3"/><text x="${left-10}" y="${assocBottom-v*assocH+4}" text-anchor="end" fill="#54675f">${(v*maxLogP).toFixed(maxLogP > 10 ? 0 : 1)}</text>`).join("");
  const legendY=height-43, legend=[0,.2,.4,.6,.8,1].map((v,i)=>`<rect x="${left+i*30}" y="${legendY}" width="30" height="12" fill="${r2Color(v)}"/><text x="${left+i*30}" y="${legendY+27}" fill="#54675f">${v}</text>`).join("");
  return `<div class="ld-metric-heading"><strong>${metricLabel} LD 热图</strong><span>${isDprime ? "|D′|（未定相基因型使用 EM 估计）" : "等位基因剂量相关平方"}</span></div><div class="ld-figure-toolbar"><span>${heatmap.plotted_count} / ${heatmap.original_count} 个连锁位点${heatmap.downsampled ? "（已按连锁强度抽样）" : ""}</span><label>PNG 倍率 <select class="ld-png-scale" data-metric="${metric}"><option>2</option><option selected>4</option><option>6</option></select></label><button class="secondary export-ld-png" data-metric="${metric}">导出高清 PNG</button><button class="secondary export-ld-pdf" data-metric="${metric}">导出矢量 PDF</button></div>
    <div class="regional-ld-figure"><svg id="regionalLdSvg-${metric}" viewBox="0 0 ${width} ${height}" role="img" aria-label="区域关联、基因结构、Lead 连锁跨度与 ${metricLabel} 三角热图"><title>CallVCF ${metricLabel} LD heatmap</title><desc>上部为区域关联信号和基因结构，中部为 Lead 连锁跨度，下部为成对 ${metricLabel}。</desc>
      <rect width="${width}" height="${height}" fill="#ffffff"/>${yTicks}${assocMarks}${leadGuides}
      <line x1="${left}" x2="${width-right}" y1="${assocBottom}" y2="${assocBottom}" stroke="#71827a"/>
      <text x="18" y="110" transform="rotate(-90 18 110)" fill="#344b43">-log10(P)</text>
      <text x="${left}" y="205" fill="#54675f">${escapeHtml(region.chrom)}:${formatNumber(region.start)}</text><text x="${width-right}" y="205" text-anchor="end" fill="#54675f">${escapeHtml(region.chrom)}:${formatNumber(region.end)}</text>
      ${geneTrack.markup}${blockMarks}<text x="${left}" y="${blockTop+62}" fill="#344b43">Lead-linked spans</text>${cells.join("")}
      <text x="${left}" y="${legendY-10}" fill="#344b43">Pairwise ${metricLabel}</text>${legend}
    </svg></div>`;
}

function renderRegionalLdFigure(data) {
  return `<div class="ld-dual-figures">${renderRegionalLdMetric(data,"r2")}${renderRegionalLdMetric(data,"dprime")}</div>`;
}

function downloadBlob(blob, name) {
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob); link.download = name; link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 1500);
}

function exportLdPng(metric="r2") {
  const svg = $(`regionalLdSvg-${metric}`);
  if (!svg) throw new Error("当前没有可导出的热图");
  const scale = Number(document.querySelector(`.ld-png-scale[data-metric="${metric}"]`)?.value || 4);
  const clone = svg.cloneNode(true), view = svg.viewBox.baseVal;
  clone.setAttribute("width", view.width * scale); clone.setAttribute("height", view.height * scale);
  const blob = new Blob([new XMLSerializer().serializeToString(clone)], {type: "image/svg+xml;charset=utf-8"});
  const image = new Image(), url = URL.createObjectURL(blob);
  image.onload = () => {
    const canvas = document.createElement("canvas"); canvas.width = view.width * scale; canvas.height = view.height * scale;
    const ctx = canvas.getContext("2d"); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, canvas.width, canvas.height); ctx.drawImage(image, 0, 0, canvas.width, canvas.height);
    URL.revokeObjectURL(url); canvas.toBlob(out => downloadBlob(out, `CallVCF_LD_${metric}_${scale}x.png`), "image/png");
  };
  image.onerror = () => { URL.revokeObjectURL(url); toast("PNG 导出失败", true); };
  image.src = url;
}

function pdfEscape(value) { return String(value).replace(/[\\()]/g, "\\$&"); }

function exportLdPdf(metric="r2") {
  const result = state.lastLeadResult, ld = result?.ld, heatmap = ld?.heatmap;
  if (!heatmap?.variants?.length) throw new Error("当前没有可导出的热图");
  const isDprime=metric==="dprime", metricLabel=isDprime?"Dprime":"R2";
  const matrix=isDprime?heatmap.matrix_dprime:(heatmap.matrix_r2||heatmap.matrix);
  const pageW = 842, pageH = 595, left = 58, right = 24, plotW = pageW-left-right;
  const heatBottom = 58, heatH = 330, n = heatmap.variants.length, cellW = plotW/n, cellH = heatH/n;
  const region = ld.search_region, regionSpan = Math.max(1, region.end-region.start), px = pos => left+(pos-region.start)/regionSpan*plotW;
  const association = (result.phenotype?.records || []).filter(r => r.chrom===region.chrom && r.pos>=region.start && r.pos<=region.end);
  const maxLogP = Math.max(1, ...association.map(r => -Math.log10(Math.max(Number(r.pvalue),1e-300))));
  const commands = ["1 1 1 rg 0 0 842 595 re f", "0.15 0.22 0.19 rg", `BT /F1 15 Tf 58 572 Td (${pdfEscape(`CallVCF regional association, genes and ${metricLabel} LD`)}) Tj ET`, `BT /F1 8 Tf 58 557 Td (${pdfEscape(`${region.chrom}:${region.start}-${region.end}; ${heatmap.plotted_count} linked variants; pairwise ${metricLabel}`)}) Tj ET`, "0.55 0.62 0.59 RG 0.6 w 58 468 m 818 468 l S"];
  association.slice(0,12000).forEach(r => {
    const y=468+(-Math.log10(Math.max(Number(r.pvalue),1e-300))/maxLogP)*68;
    commands.push(`0.35 0.45 0.41 rg ${(px(r.pos)-1).toFixed(2)} ${(y-1).toFixed(2)} 2 2 re f`);
  });
  (result.gene_track?.models||[]).slice(0,12).forEach((gene,i)=>{
    const y=432-(i%4)*10, gx1=px(gene.start), gx2=px(gene.end);
    commands.push(`0.25 0.40 0.34 RG 0.8 w ${gx1.toFixed(2)} ${y} m ${gx2.toFixed(2)} ${y} l S`);
    (gene.features||[]).filter(f=>f.type==="exon"||f.type==="cds").forEach(f=>{
      const fx1=px(f.start), fx2=px(f.end), h=f.type==="cds"?6:4;
      commands.push(`0.18 0.42 0.34 rg ${Math.min(fx1,fx2).toFixed(2)} ${(y-h/2).toFixed(2)} ${Math.max(1,Math.abs(fx2-fx1)).toFixed(2)} ${h} re f`);
    });
  });
  const leads = ld.lead_results?.length ? ld.lead_results : [{lead_record:ld.lead_record,block:ld.linkage_region}];
  leads.forEach((item,i) => {
    const y=445-(i%3)*5, block=item.block||ld.linkage_region, leadX=px(item.lead_record.pos);
    commands.push(`0.78 0.22 0.25 RG 1 w ${leadX.toFixed(2)} 420 m ${leadX.toFixed(2)} 540 l S`);
    commands.push(`0.78 0.22 0.25 RG 4 w ${px(block.start).toFixed(2)} ${y} m ${px(block.end).toFixed(2)} ${y} l S`);
  });
  commands.push("0.25 0.34 0.30 rg BT /F1 8 Tf 58 430 Td (Lead-linked spans) Tj ET");
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    const value = matrix[i]?.[j];
    const v = value == null ? 0 : Math.max(0, Math.min(1, Number(value)));
    const r = (255 - 18 * (1-v))/255, g = (246 - 202*v)/255, b = (220 - 184*v)/255;
    commands.push(`${r.toFixed(3)} ${g.toFixed(3)} ${b.toFixed(3)} rg ${(left+j*cellW).toFixed(2)} ${(heatBottom+(n-i-1)*cellH).toFixed(2)} ${(cellW+.1).toFixed(2)} ${(cellH+.1).toFixed(2)} re f`);
  }
  commands.push(`0.2 0.3 0.26 RG 0.5 w ${left} ${heatBottom} ${plotW} ${heatH} re S`);
  commands.push(`0.25 0.34 0.30 rg BT /F1 8 Tf 58 42 Td (Pairwise ${metricLabel}: 0 = pale, 1 = red) Tj ET`);
  const content = commands.join("\n");
  const objects = [
    "<< /Type /Catalog /Pages 2 0 R >>",
    "<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
    `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${pageW} ${pageH}] /Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>`,
    `<< /Length ${new TextEncoder().encode(content).length} >>\nstream\n${content}\nendstream`,
    "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
  ];
  let pdf = "%PDF-1.4\n", offsets = [0];
  objects.forEach((obj, i) => { offsets.push(new TextEncoder().encode(pdf).length); pdf += `${i+1} 0 obj\n${obj}\nendobj\n`; });
  const xref = new TextEncoder().encode(pdf).length;
  pdf += `xref\n0 ${objects.length+1}\n0000000000 65535 f \n${offsets.slice(1).map(x => String(x).padStart(10,"0")+" 00000 n ").join("\n")}\ntrailer\n<< /Size ${objects.length+1} /Root 1 0 R >>\nstartxref\n${xref}\n%%EOF`;
  downloadBlob(new Blob([pdf], {type:"application/pdf"}), `CallVCF_LD_${metric}.pdf`);
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
    const linkedRows = (ld.linked_variants || []).map(x => ({位点: `${x.chrom}:${x.pos}`, ID: x.id || ".", REF: x.ref, ALT: x.alt, 类型: x.variant_type, "最大r²": x.max_r2 ?? x.r2, 连锁Lead数: x.lead_count ?? 1, 连锁Lead: (x.linked_leads || [data.lead_locus]).join(", "), 有效样本: x.n, 距Lead_bp: x.distance ?? "—"}));
    const leadRows = (ld.lead_results || []).map(x => ({Lead: x.lead_record.key, 连锁变异数: x.linked_variant_count, 区段起点: x.block.start, 区段终点: x.block.end, 跨度_bp: x.block.span_bp}));
    html += `<section class="analysis-section"><h4>${data.ld_mode === "multi_lead_region" ? "区间内多 Lead 连锁 SNP 联集" : "指定 Lead 的连锁 SNP"}</h4><p class="sub">搜索 ${ld.search_region.chrom}:${formatNumber(ld.search_region.start)}–${formatNumber(ld.search_region.end)}；检测 ${formatNumber(ld.tested_variant_count)} 个有效二等位变异，r² ≥ ${ld.threshold}。</p>
      ${renderRegionalLdFigure(data)}
      ${leadRows.length ? `<h4 class="subheading">各 Lead 连锁跨度</h4>${genericTable(leadRows, ["Lead", "连锁变异数", "区段起点", "区段终点", "跨度_bp"])}` : ""}
      <div class="pane-intro" style="margin:12px 0 8px"><p>连锁跨度 ${formatNumber(ld.linkage_region.span_bp)} bp；r² ≥ ${ld.threshold}。</p><div><button class="secondary" id="copyLinkedBtn">复制位点</button> <button class="secondary" id="downloadLinkedBtn">下载 TSV</button></div></div>
      ${genericTable(linkedRows, ["位点", "ID", "REF", "ALT", "类型", "最大r²", "连锁Lead数", "连锁Lead", "有效样本", "距Lead_bp"])}</section>`;
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
  const ids = ["gffPath", "annotationPath", "domainPath", "phenotypePath", "ldblockshowPath", "outputDir", "ldblockshowMetric", "ldblockshowBlockType", "ldWindow", "ldThreshold", "ldMinSamples", "ldMode", "ldRegion", "leadLoci", "heatmapMaxVariants", "profilePhenotypePath", "profilePThreshold", "traitDirections"];
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

function renderToolStatus(data) {
  state.tools = data;
  const plink = data.plink, ldb = data.ldblockshow, bcf = data.bcftools;
  $("plinkStatus").textContent = plink.installed ? `已就绪 · ${plink.path}` : `未安装 · ${plink.bundled_version}`;
  $("ldblockshowStatus").textContent = ldb.installed ? (ldb.requires_wsl && !ldb.wsl_available ? `已下载 · 需先安装 WSL` : `已就绪 · ${ldb.path}`) : (ldb.requires_wsl && !ldb.wsl_available ? "未安装 · 需先启用 WSL" : "未安装");
  $("bcftoolsStatus").textContent = bcf.installed ? `已就绪 · ${bcf.path}` : (bcf.requires_wsl && !bcf.wsl_available ? "需先完成 Ubuntu/WSL 初始化" : "未安装 · 可一键安装");
  $("toolRoot").textContent = `工具目录：${data.tool_root}；PLINK GPL-3.0，LDBlockShow MIT，bcftools MIT/Expat（部分插件GPL）。`;
  document.querySelectorAll(".install-tool").forEach(button => {
    const installed = data[button.dataset.tool]?.installed;
    button.textContent = installed ? "重新安装" : "一键安装";
  });
  if (ldb.installed && !$("ldblockshowPath").value.trim()) $("ldblockshowPath").value = ldb.path;
  if (data.platform === "Windows" && ldb.installed) $("ldblockshowWsl").checked = true;
}

async function loadToolStatus() {
  try {
    const response = await fetch("/api/tools/status", {cache:"no-store"});
    const result = await response.json();
    if (!response.ok || !result.ok) throw new Error(result.error || `HTTP ${response.status}`);
    renderToolStatus(result.data);
  } catch (error) {
    $("plinkStatus").textContent = "检测失败"; $("ldblockshowStatus").textContent = "检测失败"; $("bcftoolsStatus").textContent = "检测失败";
  }
}

async function installTool(button) {
  const tool = button.dataset.tool, oldText = button.textContent;
  button.disabled = true; button.textContent = "正在下载…";
  try {
    const result = await api("/api/tools/install", {tool});
    renderToolStatus(result); toast(`${tool === "plink" ? "PLINK" : tool === "bcftools" ? "bcftools" : "LDBlockShow"} 安装完成`);
  } catch (error) { toast(error.message, true); button.textContent = oldText; }
  finally { button.disabled = false; }
}

async function runLeadAnalysis() {
  const button = $("leadRunBtn"), target = $("leadResult");
  try {
    const mode = $("ldMode").value;
    const leadText = mode === "multi_lead_region" ? $("leadLoci").value.trim() : $("leadLocus").value.trim();
    const leadLocus = leadText.split(/[\s,;]+/).find(Boolean) || lociText(false).split(/[\s,;]+/).find(Boolean);
    if (!leadLocus) throw new Error("请输入一个 Lead 位点");
    const options = { ld: $("optLd").checked, gene_track: $("optGeneTrack").checked, gene: $("optGene").checked, function: $("optFunction").checked,
      domain: $("optDomain").checked, phenotype: $("optPhenotype").checked, ldblockshow: $("optLdblockshow").checked };
    if (!Object.values(options).some(Boolean)) throw new Error("请至少勾选一项分析功能");
    if (mode === "multi_lead_region" && !$("ldRegion").value.trim()) throw new Error("多 Lead 模式请输入指定区间");
    const payload = { path: currentPath(), lead_locus: leadLocus, lead_loci: leadText || leadLocus, ld_mode: mode,
      region: $("ldRegion").value.trim(), heatmap_max_variants: $("heatmapMaxVariants").value, window_kb: $("ldWindow").value,
      r2_threshold: $("ldThreshold").value, min_samples: $("ldMinSamples").value, options,
      gff_path: $("gffPath").value.trim(), annotation_path: $("annotationPath").value.trim(),
      domain_path: $("domainPath").value.trim(), phenotype_path: $("phenotypePath").value.trim(),
      ldblockshow_path: $("ldblockshowPath").value.trim(), ldblockshow_wsl: $("ldblockshowWsl").checked,
      ldblockshow_metric: $("ldblockshowMetric").value, ldblockshow_block_type: $("ldblockshowBlockType").value,
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

function renderTrendChart(data) {
  const rows = (data.rows || []).filter(x => x.advantage_contribution != null);
  if (!rows.length) return '<div class="empty-state">没有可计算趋势的记录；请检查 β、效应等位基因和表型有利方向。</div>';
  const shown = [...rows].sort((a,b) => Math.abs(b.advantage_contribution) - Math.abs(a.advantage_contribution)).slice(0, 60);
  const width = 960, rowH = 24, top = 26, left = 210, right = 55, height = top + shown.length * rowH + 30;
  const maxAbs = Math.max(.001, ...shown.map(x => Math.abs(x.advantage_contribution))), half = (width-left-right)/2, zero = left+half;
  const marks = shown.map((x, i) => {
    const y = top+i*rowH, value = x.advantage_contribution, w = Math.abs(value)/maxAbs*half;
    const rx = value >= 0 ? zero : zero-w;
    return `<text x="${left-8}" y="${y+14}" text-anchor="end" fill="#344b43">${escapeHtml(x.sample)} · ${escapeHtml(x.trait)} · ${escapeHtml(x.lead)}</text><rect x="${rx}" y="${y+3}" width="${w}" height="14" fill="${value >= 0 ? "#23836c" : "#c94955"}"><title>贡献=${value}；P=${x.pvalue}；GT=${escapeHtml(x.gt)}</title></rect>`;
  }).join("");
  return `<div class="trend-chart"><svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Lead-SNP 表型趋势贡献图"><line x1="${zero}" x2="${zero}" y1="8" y2="${height-20}" stroke="#71827a"/><text x="${zero-half/2}" y="16" text-anchor="middle" fill="#9a3440">不利方向</text><text x="${zero+half/2}" y="16" text-anchor="middle" fill="#176b59">有利方向</text>${marks}</svg></div>`;
}

function renderSampleLeadProfile(data) {
  const cards = (data.summaries || []).map(x => `<div class="meta-card"><span>${escapeHtml(x.sample)} · 综合趋势指数</span><strong class="${x.trend_index == null ? "" : x.trend_index >= 0 ? "status-good" : "status-bad"}">${x.trend_index == null ? "未计算" : x.trend_index}</strong><small>${x.distinct_leads} 个 Lead · ${x.significant_records} 条显著 · 相对有利 ${x.favorable_records} / 不利 ${x.unfavorable_records} / 中性 ${x.neutral_records || 0}</small></div>`).join("");
  const rows = (data.rows || []).map(x => ({样本:x.sample, Lead:x.lead, 表型:x.trait, GT:x.gt, ALT剂量:x.dosage_alt, 效应等位基因:x.effect_allele || "—", 效应剂量:x.effect_copies ?? "—", 群体均值:x.cohort_mean_effect_copies ?? "—", beta:x.beta ?? "—", P值:x.pvalue, "-log10(P)":x.neg_log10_p, 显著:x.significant ? "是" : "否", 表型方向:x.trait_direction || "未指定", 等位基因效应:x.effect_direction === "favorable" ? "有利" : x.effect_direction === "unfavorable" ? "不利" : "未知", 样本相对趋势:x.trend === "favorable" ? "有利" : x.trend === "unfavorable" ? "不利" : x.trend === "neutral" ? "中性" : "方向未知", 趋势贡献:x.advantage_contribution ?? "—"}));
  const missing = data.missing_leads?.length ? `<div class="missing-box">VCF 中未找到：${data.missing_leads.map(escapeHtml).join("、")}</div>` : "";
  return `${missing}<div class="lead-summary">${cards}</div><p class="sub">${escapeHtml(data.method_note)}</p>${renderTrendChart(data)}<div class="pane-intro" style="margin:12px 0 8px"><p>匹配 ${data.association_record_count} 条 Lead×表型关联记录。</p><button class="secondary" id="downloadProfileBtn">下载明细 TSV</button></div>${genericTable(rows, ["样本","Lead","表型","GT","ALT剂量","效应等位基因","效应剂量","群体均值","beta","P值","-log10(P)","显著","表型方向","等位基因效应","样本相对趋势","趋势贡献"])}`;
}

async function runSampleLeadProfile() {
  const button = $("sampleLeadRunBtn"), target = $("sampleLeadResult");
  try {
    const payload = {path: currentPath(), samples: selectedSamples(true), lead_loci: $("profileLeadLoci").value.trim(), phenotype_path: $("profilePhenotypePath").value.trim(), significance_threshold: $("profilePThreshold").value, trait_directions: $("traitDirections").value};
    if (!payload.lead_loci) throw new Error("请输入 Lead-SNP 列表");
    setLoading(target, button, "正在联合 VCF 与 GWAS 结果生成样本画像…");
    const result = await api("/api/sample-lead-profile", payload); state.lastProfileResult = result; saveAdvancedSettings();
    clearLoading(target, button); target.innerHTML = renderSampleLeadProfile(result); toast("样本 Lead-SNP 画像已生成");
  } catch (error) { clearLoading(target, button); target.className = "result-body empty-state"; target.textContent = error.message; toast(error.message, true); }
}

function downloadProfile() {
  const rows = state.lastProfileResult?.rows || [];
  const fields = ["sample","lead","trait","gt","dosage_alt","effect_allele","effect_copies","cohort_mean_effect_copies","beta","pvalue","neg_log10_p","significant","trait_direction","effect_direction","trend","advantage_contribution"];
  const text = [fields.join("\t"), ...rows.map(row => fields.map(k => row[k] ?? "").join("\t"))].join("\n");
  downloadBlob(new Blob([text], {type:"text/tab-separated-values;charset=utf-8"}), "CallVCF_sample_lead_profile.tsv");
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

async function loadQualityCatalog() {
  try {
    const response = await fetch("/api/quality/catalog", {cache: "no-store"});
    const payload = await response.json();
    if (!response.ok || !payload.ok) throw new Error(payload.error || "无法载入质量Profile");
    state.qualityCatalog = payload.data;
    const select = $("qualityProfile");
    select.innerHTML = payload.data.profiles.map(x => `<option value="${escapeHtml(x.id)}">${escapeHtml(x.name)}</option>`).join("");
    const saved = localStorage.getItem("callvcfQualityProfile");
    select.value = payload.data.profiles.some(x => x.id === saved) ? saved : "cotton_inbred";
    applyQualityProfile();
    $("qualityOutputDir").placeholder = `默认：${payload.data.default_report_root}`;
  } catch (error) {
    toast(error.message, true);
  }
}

function applyQualityProfile() {
  const profile = (state.qualityCatalog?.profiles || []).find(x => x.id === $("qualityProfile").value);
  if (!profile) return;
  $("qualityProfileNote").textContent = profile.notes || "";
  $("qualityKingdom").value = profile.kingdom || "other";
  $("qualityPloidy").value = profile.ploidy || 2;
  $("qualityGenotypePloidy").value = String(profile.genotype_ploidy || "auto");
  $("qualitySubgenomes").value = profile.subgenomes || 1;
  $("qualityMating").value = profile.mating_system || "unknown";
  localStorage.setItem("callvcfQualityProfile", profile.id);
}

function qualityPayload() {
  const thresholds = {};
  document.querySelectorAll("[data-qc-threshold]").forEach(input => {
    if (input.value.trim() !== "") thresholds[input.dataset.qcThreshold] = input.value;
  });
  return {
    path: currentPath(),
    output_dir: $("qualityOutputDir").value.trim(),
    scan_mode: $("qualityScanMode").value,
    target_records: Number($("qualityTargetRecords").value || 200000),
    reference_path: $("qualityReferencePath").value.trim(),
    sample_meta_path: $("qualitySampleMetaPath").value.trim(),
    region_bed_path: $("qualityRegionBedPath").value.trim(),
    population_options: {
      hwe: $("populationHwe").checked,
      pca: $("populationPca").checked,
      kinship: $("populationKinship").checked,
      ld: $("populationLd").checked,
      roh: $("populationRoh").checked,
      site_missing: Number($("populationSiteMissing").value || 0.10),
      ld_max_markers: Number($("populationLdMarkers").value || 5000),
      ld_window_kb: Number($("populationLdWindow").value || 1000),
    },
    profile: {
      profile_id: $("qualityProfile").value,
      species_name: $("qualitySpecies").value.trim(),
      kingdom: $("qualityKingdom").value,
      ploidy: Number($("qualityPloidy").value),
      genotype_ploidy: $("qualityGenotypePloidy").value,
      subgenomes: Number($("qualitySubgenomes").value),
      mating_system: $("qualityMating").value,
      thresholds,
    },
  };
}

function setQualityProgress(job) {
  const wrap = $("qualityProgress");
  wrap.classList.remove("hidden");
  $("qualityProgressText").textContent = job.message || "正在运行";
  const percent = job.progress == null ? null : Math.max(0, Math.min(100, Number(job.progress)));
  $("qualityProgressValue").textContent = percent == null ? "扫描中" : `${percent.toFixed(1)}%`;
  $("qualityProgressBar").style.width = percent == null ? "35%" : `${percent}%`;
  $("qualityProgressBar").classList.toggle("indeterminate", percent == null);
  $("qualityProcessed").textContent = `已处理 ${formatNumber(job.processed_records || 0)} 条记录 · 任务 ${job.id}`;
}

function renderQualityResult(job) {
  const target = $("qualityResult");
  if (job.status === "failed" || job.status === "cancelled") {
    target.className = "result-body";
    target.innerHTML = `<div class="missing-box"><b>${job.status === "cancelled" ? "评估已取消" : "评估失败"}</b><br>${escapeHtml(job.error || job.message || "未知错误")}</div>`;
    return;
  }
  const data = job.result, summary = data.summary, scan = data.scan;
  const artifactLinks = (job.artifacts || []).map(x => {
    const open = x.name === "report.html" ? " target=\"_blank\"" : " download";
    const href = x.name === "report.html" ? `${x.url}&view=1` : x.url;
    return `<a class="artifact-link" href="${escapeHtml(href)}"${open}>${escapeHtml(x.name)} <small>${formatBytes(x.size)}</small></a>`;
  }).join("");
  const warnings = (data.warnings || []).slice(0, 200).map(x => ({
    级别: `<span class="level-${escapeHtml(x.level)}">${escapeHtml(x.level.toUpperCase())}</span>`,
    范围: x.scope, 对象: x.target, 问题: x.message, 证据: x.evidence, 建议: x.advice,
  }));
  const samples = (data.samples || []).map(x => ({
    样本: x.sample_id, Group: x.group || "—", Batch: x.batch || "—", 样本分: x.sample_score ?? "—", 状态: x.sample_status || "—", 缺失率: formatPercent(x.missing_rate), 杂合率: formatPercent(x.het_rate),
    中位DP: x.median_dp ?? "—", 中位GQ: x.median_gq ?? "—", AB异常: formatPercent(x.ab_outlier_rate),
    倍性不符: formatPercent(x.ploidy_mismatch_rate),
  }));
  const moduleLabels = {hwe:"HWE", pca:"PCA", kinship:"亲缘关系/IBD", ld:"LD衰减", roh:"ROH"};
  const populationRows = Object.entries(data.population_analysis?.modules || {}).map(([key, value]) => ({
    模块: moduleLabels[key] || key,
    状态: value.status,
    解释口径: value.interpretation,
    摘要: value.error || Object.entries(value.summary || {}).map(([k,v]) => `${k}=${v ?? "—"}`).join("；"),
    文件: (value.artifacts || []).join("、") || "—",
  }));
  target.className = "result-body";
  target.innerHTML = `<div class="quality-summary">
    <div class="meta-card"><span>质量分</span><strong>${summary.score}/100</strong></div>
    <div class="meta-card"><span>记录数</span><strong>${formatNumber(summary.record_count)}</strong></div>
    <div class="meta-card"><span>评估位点</span><strong>${formatNumber(scan.evaluated_records)}</strong></div>
    <div class="meta-card"><span>Critical</span><strong>${formatNumber(summary.critical_count)}</strong></div>
    <div class="meta-card"><span>Warning</span><strong>${formatNumber(summary.warning_count)}</strong></div>
  </div>
  <div class="notice"><b>${escapeHtml(data.profile.name)}</b> · 生物学倍性 ${data.profile.ploidy} · GT编码倍性 ${data.profile.effective_gt_ploidy ?? "未识别"} · ${scan.mode === "full" ? "完整扫描" : `智能抽样 ${formatPercent(scan.sampling_fraction)}`} · ${escapeHtml(scan.method_note)}</div>
  <div class="quality-section"><h4>HWE、PCA、亲缘关系、LD与ROH</h4><p class="sub">HWE按物种Profile解释；棉花/自交/多倍体只做描述。PCA/亲缘使用LD剪枝面板，LD衰减使用未剪枝QC面板；亲缘估计为PLINK 1.9 IBD/PI_HAT。</p>${genericTable(populationRows, ["模块","状态","解释口径","摘要","文件"])}</div>
  <div class="quality-section"><h4>报告下载</h4><div class="artifact-list">${artifactLinks}</div><small>HTML可离线交互和打印为PDF；ZIP包含HTML、JSON、TSV与运行审计文件。</small></div>
  <div class="quality-section"><h4>告警明细（页面最多显示200条，完整内容见 warnings.tsv）</h4>${genericTable(warnings, ["级别","范围","对象","问题","证据","建议"])}</div>
  <div class="quality-section"><h4>样本指标</h4>${genericTable(samples, ["样本","Group","Batch","样本分","状态","缺失率","杂合率","中位DP","中位GQ","AB异常","倍性不符"])}</div>`;
}

async function pollQuality(runId) {
  clearTimeout(state.qualityPoll);
  try {
    const job = await api("/api/quality/status", {run_id: runId});
    state.qualityRun = job;
    setQualityProgress(job);
    if (["complete", "failed", "cancelled"].includes(job.status)) {
      $("qualityRunBtn").disabled = false;
      $("qualityRunBtn").textContent = "重新评估";
      $("qualityCancelBtn").classList.add("hidden");
      renderQualityResult(job);
      if (job.status === "complete") toast("VCF质量报告已生成");
      return;
    }
    state.qualityPoll = setTimeout(() => pollQuality(runId), 900);
  } catch (error) {
    $("qualityRunBtn").disabled = false;
    $("qualityCancelBtn").classList.add("hidden");
    toast(error.message, true);
  }
}

async function runQualityAssessment() {
  try {
    const button = $("qualityRunBtn");
    button.disabled = true; button.textContent = "正在启动…";
    $("qualityCancelBtn").classList.remove("hidden");
    $("qualityResult").className = "result-body loading";
    $("qualityResult").textContent = "正在创建质量评估任务…";
    const job = await api("/api/quality/start", qualityPayload());
    state.qualityRun = job;
    button.textContent = "评估中…";
    setQualityProgress(job);
    pollQuality(job.id);
  } catch (error) {
    $("qualityRunBtn").disabled = false; $("qualityRunBtn").textContent = "开始评估";
    $("qualityCancelBtn").classList.add("hidden");
    $("qualityResult").className = "result-body empty-state";
    $("qualityResult").textContent = "尚未运行质量评估。";
    toast(error.message, true);
  }
}

async function cancelQualityAssessment() {
  if (!state.qualityRun?.id) return;
  try { await api("/api/quality/cancel", {run_id: state.qualityRun.id}); toast("正在取消质量评估"); }
  catch (error) { toast(error.message, true); }
}

function showDangerRepairPolicy(plan = null) {
  state.repairPlan = plan;
  $("dangerConfirmText").value = "";
  if (plan) {
    $("dangerRepairTitle").textContent = "危险修复必须二次确认";
    $("dangerPlanSummary").innerHTML = `<b>${escapeHtml(plan.action_name)}</b><br>输入：${escapeHtml(plan.source)}<br>输出：${escapeHtml(plan.output || "—")}<br>命令：${escapeHtml((plan.command_preview || []).join(" "))}`;
    $("dangerConfirmLabel").textContent = `请输入：${plan.confirmation_phrase}`;
    $("dangerConfirmBtn").textContent = "确认并执行";
    $("dangerConfirmBtn").disabled = true;
  } else {
    $("dangerRepairTitle").textContent = "自动修复安全规则";
    $("dangerPlanSummary").innerHTML = "<b>永不支持：</b>覆盖原VCF或已有目标、删除文件、静默改写GT。标准化和过滤必须逐次生成计划并输入本次专属短语。";
    $("dangerConfirmLabel").textContent = "此窗口仅说明规则，无修复计划";
    $("dangerConfirmBtn").textContent = "关闭";
    $("dangerConfirmBtn").disabled = false;
  }
  $("dangerRepairModal").classList.remove("hidden");
  if (plan) $("dangerConfirmText").focus();
}

function closeDangerRepairPolicy() { $("dangerRepairModal").classList.add("hidden"); }

async function loadRepairCatalog() {
  try {
    const response = await fetch("/api/repair/catalog", {cache:"no-store"});
    const payload = await response.json();
    if (!response.ok || !payload.ok) throw new Error(payload.error || "无法检测修复后端");
    state.repairCatalog = payload.data;
    $("repairAvailability").textContent = payload.data.available
      ? `bcftools已就绪：${payload.data.backend}。所有副本操作均写入新文件。`
      : "未检测到bcftools：修复执行器不可用；质量评估与PLINK群体分析不受影响。";
    $("repairPlanBtn").disabled = !payload.data.available;
  } catch (error) { $("repairAvailability").textContent = error.message; $("repairPlanBtn").disabled = true; }
}

function updateRepairFields() {
  const action = $("repairAction").value;
  $("repairExecutorCard").querySelector(".repair-output-field").classList.toggle("hidden", action === "index");
  $("repairExecutorCard").querySelector(".repair-reference-field").classList.toggle("hidden", action !== "normalize_copy");
  $("repairExecutorCard").querySelector(".repair-expression-field").classList.toggle("hidden", !["filter_copy", "mask_genotypes_copy"].includes(action));
  $("repairExecutorCard").querySelector(".repair-samples-field").classList.toggle("hidden", action !== "subset_samples_copy");
  if (action === "mask_genotypes_copy") {
    $("repairExpressionLabel").textContent = "要掩蔽为缺失的GT条件";
    $("repairExpression").placeholder = "例如：FMT/DP<5 || FMT/GQ<20";
  } else {
    $("repairExpressionLabel").textContent = "要保留位点的bcftools表达式";
    $("repairExpression").placeholder = "例如：QUAL>=30 && F_MISSING<0.1";
  }
}

function renderRepairPlan(plan) {
  $("repairResult").innerHTML = `<b>${escapeHtml(plan.action_name)}</b> · 风险：${escapeHtml(plan.risk)}<br>输入：${escapeHtml(plan.source)}<br>输出：${escapeHtml(plan.output || "仅生成索引旁文件")}<br><code>${escapeHtml((plan.command_preview || []).join(" "))}</code>`;
}

async function executeRepair(plan, confirmation = "") {
  closeDangerRepairPolicy();
  const job = await api("/api/repair/execute", {plan_id: plan.id, confirmation});
  $("repairPlanBtn").disabled = true; $("repairCancelBtn").classList.remove("hidden");
  pollRepair(job.id);
}

async function createRepairPlan() {
  try {
    const plan = await api("/api/repair/plan", {
      action: $("repairAction").value, path: currentPath(), output_path: $("repairOutputPath").value.trim(),
      reference_path: $("repairReferencePath").value.trim(), expression: $("repairExpression").value.trim(), samples: $("repairSamples").value.trim(),
    });
    state.repairPlan = plan; renderRepairPlan(plan);
    if (plan.risk === "dangerous") showDangerRepairPolicy(plan);
    else if (window.confirm(`核对完成后执行“${plan.action_name}”？\n原VCF不会被覆盖。`)) await executeRepair(plan);
  } catch (error) { toast(error.message, true); }
}

async function pollRepair(planId) {
  clearTimeout(state.repairPoll);
  try {
    const job = await api("/api/repair/status", {plan_id: planId});
    $("repairResult").innerHTML = `<b>${escapeHtml(job.message)}</b> · ${job.progress ?? 0}%${job.error ? `<br>${escapeHtml(job.error)}` : ""}${job.result ? `<br>输出：${escapeHtml(job.result.output_path || job.result.index_path)}<br>审计日志：${escapeHtml(job.result.audit_log)}${job.result.before_stats ? `<br>修复前后统计：${escapeHtml(job.result.before_stats)} · ${escapeHtml(job.result.after_stats)}<br>清单：${escapeHtml(job.result.manifest)}` : ""}` : ""}`;
    if (["complete","failed","cancelled"].includes(job.status)) {
      $("repairPlanBtn").disabled = !state.repairCatalog?.available; $("repairCancelBtn").classList.add("hidden");
      toast(job.status === "complete" ? "安全修复已完成" : job.message, job.status !== "complete");
      return;
    }
    state.repairPoll = setTimeout(() => pollRepair(planId), 800);
  } catch (error) { toast(error.message, true); }
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
  $("sampleLeadRunBtn").addEventListener("click", runSampleLeadProfile);
  $("qualityRunBtn").addEventListener("click", runQualityAssessment);
  $("qualityCancelBtn").addEventListener("click", cancelQualityAssessment);
  $("qualityProfile").addEventListener("change", applyQualityProfile);
  $("showRepairPolicyBtn").addEventListener("click", () => showDangerRepairPolicy());
  $("repairAction").addEventListener("change", updateRepairFields);
  $("repairPlanBtn").addEventListener("click", createRepairPlan);
  $("repairCancelBtn").addEventListener("click", async () => {
    if (!state.repairPlan?.id) return;
    try { await api("/api/repair/cancel", {plan_id: state.repairPlan.id}); } catch (error) { toast(error.message, true); }
  });
  $("dangerCancelBtn").addEventListener("click", closeDangerRepairPolicy);
  $("dangerConfirmText").addEventListener("input", e => {
    $("dangerConfirmBtn").disabled = !!state.repairPlan && e.target.value.trim() !== state.repairPlan.confirmation_phrase;
  });
  $("dangerConfirmBtn").addEventListener("click", async () => {
    if (!state.repairPlan) return closeDangerRepairPolicy();
    try { await executeRepair(state.repairPlan, $("dangerConfirmText").value.trim()); }
    catch (error) { toast(error.message, true); }
  });
  $("dangerRepairModal").addEventListener("click", e => { if (e.target === $("dangerRepairModal")) closeDangerRepairPolicy(); });
  $("ldMode").addEventListener("change", () => {
    const multi = $("ldMode").value === "multi_lead_region";
    document.querySelectorAll(".multi-lead-field").forEach(x => x.classList.toggle("hidden", !multi));
    $("leadLocus").closest("label").classList.toggle("hidden", multi);
    $("ldWindow").closest("label").classList.toggle("hidden", multi);
  });
  document.querySelectorAll(".browse-resource").forEach(btn => btn.addEventListener("click", () => browseResource(btn)));
  document.querySelectorAll(".install-tool").forEach(btn => btn.addEventListener("click", () => installTool(btn)));
  $("leadResult").addEventListener("click", e => {
    if (e.target.closest("#copyLinkedBtn")) copyLinked().catch(err => toast(err.message, true));
    if (e.target.closest("#downloadLinkedBtn")) downloadLinked();
    const pngButton=e.target.closest(".export-ld-png"), pdfButton=e.target.closest(".export-ld-pdf");
    if (pngButton) { try { exportLdPng(pngButton.dataset.metric); } catch (err) { toast(err.message, true); } }
    if (pdfButton) { try { exportLdPdf(pdfButton.dataset.metric); } catch (err) { toast(err.message, true); } }
  });
  $("sampleLeadResult").addEventListener("click", e => { if (e.target.closest("#downloadProfileBtn")) downloadProfile(); });
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
  loadToolStatus();
  loadQualityCatalog();
  loadRepairCatalog();
  updateRepairFields();
  $("ldMode").dispatchEvent(new Event("change"));
  checkHealth();
});
