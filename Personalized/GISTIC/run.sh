#01 make in.file
sh run_GISTIC_01.sh &&
#02 run cosmic driver
sh run_COSMIC_driver.sh &&

#画图准备
sh a-FigureFile_OCCC.sh &&

#huatu 
sh run_gisticplot.sh &&

#combin
sh combin.sh
