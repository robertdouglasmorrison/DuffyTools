# violinPlot.R -- wrapper around the ggplot geom_violin plot function

violinPlot <- function( df, aes, ..., horiz=FALSE, boxwid=NULL, facet=NULL,
			main="", xlab="X", ylab="Y", log="") {

	require( ggplot2)

	p <- ggplot( df, aes) + geom_violin( ...)

	if (horiz) p <- p + coord_flip()

	if ( ! is.null( facet)) {
		p <- p + facet_wrap( paste( "~", as.character(facet)), scales="free_x")
	}

	# add boxplot?
	if ( ! is.null(boxwid)) {
		p <- p + geom_boxplot( width=boxwid)
	}

	# any X,Y log scaling?
	if ( regexpr( 'Y', toupper(log)) > 0) p <- p + scale_y_log10()
	if ( regexpr( 'X', toupper(log)) > 0) p <- p + scale_x_log10()

	p <- p + labs( title=main, x=xlab, y=ylab)

	p
}
