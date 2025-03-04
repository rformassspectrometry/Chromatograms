#' Run devtools::test_active_file(file = "tests/testthat/test_Chromatograms-plot.R")
#' if make changes to this file.

test_that("plotChromatograms works", {
    vdiffr::expect_doppelganger(
        "plotChromatograms-single",
        function() plotChromatograms(c_full[1])
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-multiple",
        function() plotChromatograms(c_full)
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-color",
        function() plotChromatograms(c_full[1:2], col = c("green", "blue"))
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-one-color",
        function() plotChromatograms(c_full[1:2], col = c("green"))
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-toomany-color",
        function() plotChromatograms(c_full[1:3], col = c("green", "blue"))
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-toomany-main",
        function() plotChromatograms(c_full[1:3], main = c("test1", "test2"))
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-one-title",
        function() plotChromatograms(c_full[1:2], main = "Test Title")
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-asp05",
        function() plotChromatograms(c_full, asp = 1/2)
    )

    vdiffr::expect_doppelganger(
        "plotChromatograms-asp2",
        function() plotChromatograms(c_full, asp = 2)
    )
})

test_that("plotChromatogramsOverlay works", {
    vdiffr::expect_doppelganger(
        "plotChromatogramsOverlay-basic",
        function() plotChromatogramsOverlay(c_full, col = c("red", "blue"))
    )

    vdiffr::expect_doppelganger(
        "plotChromatogramsOverlay-xlim",
        function() plotChromatogramsOverlay(c_full, xlim = c(10, 50))
    )

    vdiffr::expect_doppelganger(
        "plotChromatogramsOverlay-no-axes",
        function() plotChromatogramsOverlay(c_full, axes = FALSE)
    )

    vdiffr::expect_doppelganger(
        "plotChromatogramsOverlay-main-title",
        function() plotChromatogramsOverlay(c_full, main = "Overlay Test")
    )

        vdiffr::expect_doppelganger(
            "plotChromatogramsOverlay-one-sample",
            function() plotChromatogramsOverlay(c_full[1], col = c("red"))
        )
})

test_that(".plot_single_chromatogram works", {
    vdiffr::expect_doppelganger(
        "plot_single_chromatogram-basic",
        function() .plot_single_chromatogram(c_full[1])
    )

    vdiffr::expect_doppelganger(
        "plot_single_chromatogram-xlim",
        function() .plot_single_chromatogram(c_full[1], xlim = c(10, 50),
                                             ylim = c(0, 3000))
    )

    vdiffr::expect_doppelganger(
        "plot_single_chromatogram-color",
        function() .plot_single_chromatogram(c_full[1], col = "purple")
    )

    vdiffr::expect_doppelganger(
        "plot_single_chromatogram-infinite-xlim",
        function() .plot_single_chromatogram(c_full[1], xlim = c(-Inf, Inf))
    )

    vdiffr::expect_doppelganger(
        "plot_single_chromatogram-infinite-ylim",
        function() .plot_single_chromatogram(c_full[1], ylim = c(-Inf, Inf))
    )
})
