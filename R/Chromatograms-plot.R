#' @title Plot chromatograms
#'
#' @name plotChromatograms
#'
#' @description
#'
#' [Chromatograms()] can be plotted with the following functions:
#'
#' The `plotChromatograms()`: plots each chromatogram in its separate plot by
#' splitting the plot area into as many panels as there are spectra.
#'
#' @param x A [Chromatograms] object.
#'
#' @param xlab `character(1)` with the label for the x-axis (by default
#'        `xlab = "rtime (s)"`).
#'
#' @param ylab `character(1)` with the label for the y-axis (by default
#'        `ylab = "intensity"`).
#'
#' @param type `character(1)` specifying the type of plot. See [plot.default()]
#'        for details. Defaults to `type = "l"` which draws each peak as a line.
#'
#' @param xlim `numeric(2)` defining the x-axis limits. The range of m/z values
#'        are used by default.
#'
#' @param ylim `numeric(2)` defining the y-axis limits. The range of intensity
#'        values are used by default.
#'
#' @param main `character(1)` with the title for the plot. By default the
#'        spectrum's MS level and retention time (in seconds) is used.
#'
#' @param col color to be used to draw the peaks. Should be either of length 1,
#'        or equal to the number of chromatograms (to plot each chromatograms
#'        in a different color) or be a `list` with colors for each individual
#'        peak in each spectrum.
#'
#' @param axes `logical(1)` whether (x and y) axes should be drawn.
#'
#' @param asp `numeric(1)` the aspect ratio of the plot, i.e. the ratio of
#'        the y-axis to the x-axis. Defaults to 1.
#'
#' @param pch `numeric(1)` specifying the symbol to be used for the peaks.
#'        Defaults to 20, a filled circle. See [points()] for details.
#'
#' @param cex `numeric(1)` specifying the character expansion factor for the
#'        peaks. Defaults to 0.6.
#'
#' @param lwd `numeric(1)` specifying the line width for the peaks. Defaults to
#'        1.5.
#'
#' @param frame.plot `logical(1)` whether a box should be drawn around the
#'        plotting area.
#'
#' @param ... Additional arguments to be passed to [plot.default()].
#'
#' @return These functions create a plot.
#'
#' @author Philippine Louail, Johannes Rainer
#'
#' @examples
#'
#' ## Create a Chromatograms object
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     chromIndex = c(1L, 2L, 3L)
#' )
#' pdata <- list(
#'     data.frame(
#'         rtime = c(12.4, 12.8, 13.2, 14.6),
#'         intensity = c(123.3, 153.6, 2354.3, 243.4)
#'     ),
#'     data.frame(
#'         rtime = c(45.1, 46.2),
#'         intensity = c(100, 80.1)
#'     ),
#'     data.frame(
#'         rtime = c(12.4, 12.8, 13.2, 14.6),
#'         intensity = c(123.3, 153.6, 2354.3, 243.4)
#'     )
#' )
#' chr <- backendInitialize(ChromBackendMemory(),
#'     chromData = cdata,
#'     peaksData = pdata
#' ) |> Chromatograms()
#'
#' ## Plot one chromatogram
#' plotChromatograms(chr[1])
#'
#' ## Plot the full Chromatograms object
#' plotChromatograms(chr)
#'
#' ## Define a color for each peak in each chromatogram
#' plotChromatograms(chr[1:2], col = c("green", "blue"))
#'
#' ## Overlay all chromatograms
#' plotChromatogramsOverlay(chr[1:2], col = c("green", "blue"))
#'
NULL


#' @rdname plotChromatograms
#' @importFrom graphics par
#' @importFrom grDevices n2mfrow
#' @export
plotChromatograms <- function(x, xlab = "rtime (s)", ylab = "intensity",
    type = "o",
    pch = 20, cex = 0.6, lwd = 1.5,
    xlim = numeric(), ylim = numeric(),
    main = character(), col = "#00000080",
    asp = 1, ...) {
    if (!length(main)) {
        main <- paste0("m/z: ", round(mz(x), 1))
    } # maybe range is better
    nsp <- length(x)
    if (nsp == 1) {
        col <- list(col)
    }
    if (length(col) != nsp) {
        col <- rep(col[1], nsp)
    }
    if (length(main) != nsp) {
        main <- rep(main[1], nsp)
    }
    if (nsp > 1) {
        par(mfrow = n2mfrow(nsp, asp = asp))
    }
    for (i in seq_len(nsp)) {
        .plot_single_chromatogram(x[i],
            xlab = xlab, ylab = ylab, type = type,
            xlim = xlim, ylim = ylim, main = main[i],
            col = col[[i]], pch = pch, cex = cex,
            lwd = lwd, ...
        )
    }
}

#' @rdname plotChromatograms
#' @export
plotChromatogramsOverlay <- function(x,
    xlab = "rtime (s)", ylab = "intensity",
    type = "o",
    pch = 20, cex = 0.6, lwd = 1.5,
    xlim = numeric(),
    ylim = numeric(),
    main = paste(length(x), "chromatograms"),
    col = "#00000080",
    axes = TRUE, frame.plot = axes, ...) {
    nsp <- length(x)
    if (nsp == 1) {
        col <- list(col)
    }
    if (length(col) != nsp) {
        col <- rep(col[1], nsp)
    }
    if (!length(xlim)) {
        xlim <- range(unlist(rtime(x)), na.rm = TRUE)
    }
    if (!length(ylim)) {
        ylim <- c(0, max(unlist(intensity(x)), na.rm = TRUE))
    }
    dev.hold()
    on.exit(dev.flush())
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    if (axes) {
        axis(side = 1, ...)
        axis(side = 2, ...)
    }
    if (frame.plot) {
        box(...)
    }
    title(main = main, xlab = xlab, ylab = ylab, ...)
    for (i in seq_len(nsp)) {
        .plot_single_chromatogram(x[i],
            add = TRUE, type = type,
            col = col[[i]], pch = pch, cex = cex,
            lwd = lwd,
            ...
        )
    }
}
