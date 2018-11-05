# x11Tools.R --  routines that deal with X11 graphics


`checkX11` <- function( bg='white', type="dbcairo", width=10, height=7, ...) {

	curDevName <- names( dev.cur())

	if ( regexpr( "X11", curDevName) < 1) {
		# current device is not X11
		X11( type=type, bg=bg, width=width, height=height, ...)
		curDevName <- names( dev.cur())
	}
	curDevName
}
