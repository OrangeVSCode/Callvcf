EMMAX bundled runtime
=====================

Upstream: https://csg.sph.umich.edu/kang/emmax/download/index.html
Distribution: emmax-intel-binary-20120210.tar.gz
Archive SHA256: E2A582851BA1BE908757D4EF436E98AD76664A0C55E00D13E55FA35FE2BA54DD
Platform: Ubuntu x86_64 (invoked through WSL on Windows)
License: MIT

Files:
- emmax-intel64 SHA256 AE2A0FFA285B846B785EFE68F0CEFEB9BCD4B2A288E53E8368FD6EF01DD8A76B
- emmax-kin-intel64 SHA256 CCE230D866AB2A87DA224B84D742493584D91FF65614D38FAE9C99D31C6B8F9F

GPA-Accelerator deploys these binaries to its managed tool directory. Windows
still requires an initialized WSL Linux distribution because upstream does not
publish a native Windows executable.
