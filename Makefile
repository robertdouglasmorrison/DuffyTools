#
# DuffyTools  ver 1.8.1

VER = 1.8.1

DuffyTools.tar.gz : 

	rm  -f ./DuffyTools*.tar.gz
	${R_PROGRAM}  CMD SHLIB src/*.c
	${R_PROGRAM}  CMD build  --force .
	${R_PROGRAM}  CMD INSTALL ./DuffyTools_${VER}.tar.gz

