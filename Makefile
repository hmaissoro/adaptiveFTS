# Makefile for adaptiveFTS — common development and CRAN-submission tasks.
#
# Usage: `make <target>` (requires R, Rscript and the C++ toolchain on PATH).
# Override the interpreters if they are not on PATH, e.g.
#   make check R="/c/Program Files/R/R-4.2.1/bin/x64/R.exe" \
#              RSCRIPT="/c/Program Files/R/R-4.2.1/bin/x64/Rscript.exe"

PKG     := $(shell sed -n 's/^Package: //p' DESCRIPTION)
VERSION := $(shell sed -n 's/^Version: //p' DESCRIPTION)
TARBALL := $(PKG)_$(VERSION).tar.gz
R       ?= R
RSCRIPT ?= Rscript

.PHONY: help deps document test build check check-cran install coverage \
        lint vignettes pkgdown clean cran submit

help:
	@echo "adaptiveFTS $(VERSION) — make targets:"
	@echo "  deps        install development dependencies"
	@echo "  document    regenerate Rd files and NAMESPACE with roxygen2"
	@echo "  test        run the testthat suite (devtools::test)"
	@echo "  build       build the source tarball ($(TARBALL))"
	@echo "  install     install the package from source"
	@echo "  check       R CMD check --as-cran on the built tarball"
	@echo "  cran        document + build + check --as-cran (pre-submission gate)"
	@echo "  coverage    test-coverage report (covr)"
	@echo "  lint        lintr::lint_package() (informational)"
	@echo "  vignettes   build the vignettes"
	@echo "  pkgdown     build the documentation website"
	@echo "  submit      submit to CRAN (devtools::submit_cran; interactive)"
	@echo "  clean       remove build artefacts"

deps:
	$(RSCRIPT) -e 'pkgs <- c("devtools","roxygen2","testthat","covr","lintr","pkgdown","rcmdcheck","rhub"); ip <- rownames(installed.packages()); for (p in pkgs) if (!(p %in% ip)) install.packages(p)'

document:
	$(RSCRIPT) -e 'roxygen2::roxygenise()'

test:
	$(RSCRIPT) -e 'devtools::test()'

build: document
	$(R) CMD build .

install:
	$(R) CMD INSTALL .

check: build
	$(R) CMD check --as-cran "$(TARBALL)"

# Full local pre-submission gate: regenerate docs, build, and run --as-cran.
cran: document
	$(R) CMD build .
	$(R) CMD check --as-cran "$(TARBALL)"

coverage:
	$(RSCRIPT) -e 'covr::report(covr::package_coverage())'

lint:
	$(RSCRIPT) -e 'print(lintr::lint_package())'

vignettes:
	$(RSCRIPT) -e 'devtools::build_vignettes()'

pkgdown:
	$(RSCRIPT) -e 'pkgdown::build_site()'

submit:
	$(RSCRIPT) -e 'devtools::submit_cran()'

clean:
	rm -rf $(PKG).Rcheck $(PKG)_*.tar.gz docs
	rm -f src/*.o src/*.so src/*.dll
