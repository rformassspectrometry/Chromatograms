# `testthat` unit test suites for `ChromBackend`

This folder contains unit test suites that can be called by other packages
defining `ChromBackend` instances to ensure that they are compliant with the
`ChromBackend` definition and that they provide all necessary
methods/functionality.

Different test suites are available within sub-folders of this directory, each
of them having one or multiple *.R* files defining the unit tests. A description
of the unit tests and which variables need to be defined by the calling package
is provided within these.

The unit tests from one test suite can be called from within another package (in
their `testthat.R` file) like shown in the example below that will call all unit
tests from the *test_ChromBackend* test suite/directory. This tests perform all
tests on a backend class assigned to a variable `be`, thus developers should
ensure to have create an instance of their `ChromBackend` implementation and assign
it to a variable `be` before the code lines below.

```
test_suite <- system.file("test_backends", "test_ChromBackend",
    package = "Chromatograms")
test_dir(test_suite, stop_on_failure = TRUE)
```

Note that setting `stop_on_failure = TRUE` is required since the tests would
otherwise silently fail.  Alternatively, it is also possible to run single test
files, e.g. the test file that checks that all spectra variables can be
correctly accessed or set:

```
test_file(file.path(test_suite, "test_peaksData.R"),
    stop_on_failure = TRUE)
```

## Authors/contributors

- Philippine Louail (@philouail)
- Johannes Rainer (@jorainer)
