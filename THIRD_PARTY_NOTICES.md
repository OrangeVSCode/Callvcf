# Optional third-party tools

GPA-Accelerator can download and execute the following independent command-line tools at the user's request. They are not part of GPA-Accelerator's source code and are stored separately in the user's local application-data directory.

## PLINK 1.9

- Project: https://www.cog-genomics.org/plink/1.9/
- Version selected by GPA-Accelerator: beta 7.11, 19 August 2025
- License: GNU General Public License version 3
- GPA-Accelerator preserves the license file included in the official PLINK archive.

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

The tools remain separate programs communicating with GPA-Accelerator through files and command-line arguments. Their authors do not endorse GPA-Accelerator, and GPA-Accelerator does not provide a warranty for third-party software.
