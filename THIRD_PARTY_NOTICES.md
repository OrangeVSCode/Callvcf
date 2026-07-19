# Optional third-party tools

CallVCF can download and execute the following independent command-line tools at the user's request. They are not part of CallVCF's source code and are stored separately in the user's local application-data directory.

## PLINK 1.9

- Project: https://www.cog-genomics.org/plink/1.9/
- Version selected by CallVCF: beta 7.11, 19 August 2025
- License: GNU General Public License version 3
- CallVCF preserves the license file included in the official PLINK archive.

## LDBlockShow

- Project: https://github.com/hewm2008/LDBlockShow
- License: MIT
- Platform note: the upstream project supports Linux, Unix and macOS. CallVCF invokes it through WSL on Windows.
- CallVCF preserves the complete upstream archive, including its license and bundled helper files.

The tools remain separate programs communicating with CallVCF through files and command-line arguments. Their authors do not endorse CallVCF, and CallVCF does not provide a warranty for third-party software.
