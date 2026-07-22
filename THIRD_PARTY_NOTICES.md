# Third-party tools

GPA-Accelerator uses the following independent command-line tools. PLINK, LDBlockShow and BCFtools are installed on demand. The official EMMAX binary distribution is bundled with the complete Windows package and is deployed to the user's local application-data directory without a network download.

## PLINK 1.9

- Project: https://www.cog-genomics.org/plink/1.9/
- Version selected by GPA-Accelerator: beta 7.11, 19 August 2025
- License: GNU General Public License version 3
- GPA-Accelerator preserves the license file included in the official PLINK archive.

## PLINK 2

- Project: https://www.cog-genomics.org/plink/2.0/
- Version selected by GPA-Accelerator: alpha 7.1, 4 May 2026
- License: GNU General Public License version 3
- Purpose in GPA-Accelerator: KING-robust kinship coefficient and `.kin0` output for the sample-similarity workflow.
- GPA-Accelerator downloads the official platform archive on demand and preserves its included license file.

## LDBlockShow

- Project: https://github.com/hewm2008/LDBlockShow
- License: MIT
- Platform note: the upstream project supports Linux, Unix and macOS. GPA-Accelerator invokes it through WSL on Windows.
- GPA-Accelerator preserves the complete upstream archive, including its license and bundled helper files.

## BCFtools / HTSlib

- Project: https://github.com/samtools/bcftools
- GPA-Accelerator can install the Ubuntu distribution package into the user's WSL environment; it is not silently downloaded during application startup.
- BCFtools and HTSlib use the MIT/Expat license; optional GPL plugins/features remain subject to their applicable licenses.
- GPA-Accelerator records the managed backend and only reports success after bcftools exits successfully and the new output plus index pass validation.

## EMMAX

- Project and official download: https://csg.sph.umich.edu/kang/emmax/download/index.html
- Bundled distribution: `emmax-intel-binary-20120210.tar.gz`; program build `emmax-intel64-20120205`
- Upstream archive SHA-256: `E2A582851BA1BE908757D4EF436E98AD76664A0C55E00D13E55FA35FE2BA54DD`
- License: MIT; a copy is preserved in `vendor/emmax/LICENSE.txt` and in the managed deployment directory.
- Platform note: the official bundled executables are Ubuntu x86-64 binaries. On Windows, GPA-Accelerator invokes them through an initialized Ubuntu/WSL environment. The application bundle does not silently install or modify the Windows optional WSL feature.

The tools remain separate programs communicating with GPA-Accelerator through files and command-line arguments. Their authors do not endorse GPA-Accelerator, and GPA-Accelerator does not provide a warranty for third-party software.
