;; To use this file to build HEAD of smoothxg:
;;
;;   guix build -f guix.scm
;;
;; (make sure you are running a recent guix)
;;
;; To do a cross compilation build for ARM64
;;
;;   guix build -f guix.scm --target=aarch64-linux
;;
;; To get a development container (inside emacs shell will work)
;;
;;   guix shell -C -D -f guix.scm -- bash --init-file <(echo "ln -s /bin/sh /bin/bash")
;;
;; and build
;;
;;   find -name CMakeCache.txt|xargs rm -v
;;   cd build
;;   cmake -DCMAKE_BUILD_TYPE=Debug ..
;;   cmake --build . --verbose -- -j 14 && ctest . --verbose
;;
;; For the tests you may need /usr/bin/env. In a container create it with
;;
;;   mkdir -p /usr/bin ; ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env
;;

(use-modules
  (ice-9 popen)
  (ice-9 rdelim)
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix download)
  (guix git-download)
  (guix build-system cmake)
  (guix utils)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages commencement) ; gcc-toolchain
  (gnu packages curl)
  (gnu packages datastructures)
  (gnu packages gdb)
  (gnu packages gcc)
  (gnu packages jemalloc)
  (gnu packages libffi)
  (gnu packages mpi)
  (gnu packages python)
  (gnu packages python-xyz)
  (gnu packages pkg-config)
  (gnu packages tls)
  (gnu packages version-control)
)

(define %source-dir (dirname (current-filename)))

(define %git-commit
  (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public smoothxg-git
  (package
    (name "smoothxg-git")
    (version (git-version "0.3.3" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (inputs
     `(
       ("coreutils" ,coreutils)
       ("pybind11" ,pybind11) ;; see libstd++ note in remarks above
       ("jemalloc" ,jemalloc)
       ("gcc" ,gcc-11)
       ("gcc-lib" ,gcc-11 "lib")
       ("gcc-toolchain" ,gcc-toolchain)
       ("gdb" ,gdb)
       ("git" ,git) ; pulls in perl which does not do RISV-V cross builds yet
       ("openmpi" ,openmpi)
       ("python" ,python)
       ("sdsl-lite" ,sdsl-lite)
       ("libdivsufsort" ,libdivsufsort)
       ("zlib" ,zlib)
       ("zstd-lib" ,zstd "lib")
       ))
    (native-inputs
     `(("pkg-config" ,pkg-config)
       ))
    (arguments
      `(#:phases
        (modify-phases
         %standard-phases
         ;; This stashes our build version in the executable
         (add-after 'unpack 'set-version
           (lambda _
             (mkdir-p "include")
             (with-output-to-file "include/smoothxg_git_version.hpp"
               (lambda ()
                 (format #t "#define ODGI_GIT_VERSION \"~a\"~%" version)))
             #t))
         ;; (delete 'check)
         )
        ;; #:make-flags (list ,(string-append "CC=" (cc-for-target)))))
        ))
     (synopsis "odgi pangenome optimized dynamic sequence graph implementation")
     (description
"odgi pangenome graph tooling provides an efficient, succinct dynamic
DNA sequence graph model, as well as a host of algorithms that allow
the use of such graphs.")
     (home-page "https://github.com/vgteam/odgi")
     (license license:expat)))

smoothxg-git
