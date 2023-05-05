#
# DuffyTools  ver 1.9.2

VER = 1.9.2

DuffyTools.tar.gz : 

	rm  -f ./DuffyTools*.tar.gz
	${R_PROGRAM}  CMD SHLIB src/*.c
	${R_PROGRAM}  CMD build  --force .
	${R_PROGRAM}  CMD INSTALL ./DuffyTools_${VER}.tar.gz

