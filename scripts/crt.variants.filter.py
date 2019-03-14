#!/bin/env python

# reports variants only in muscular genes using csv report on rare variants

import csv
import sys
import re

original_report = sys.argv[1]
panel_report = original_report.replace("csv","panel.csv")


with open(original_report,'rb') as f_original_report:
    header_line = f_original_report.readline().strip().replace('"','')
    header = header_line.split(',')
f_original_report.close()

with open(original_report,'rb') as f_original_report:
	reader = csv.DictReader(f_original_report)
	with open(panel_report,'w') as f_panel_report:
	    muscle_genes=['ACTA1','ACVR1','AGRN','ALG13','ALG14','ALG2','ANO5','ATP2A1','B3GALNT2','B3GNT1','BAG3','BIN1','CACNA1A','CACNA1S','CAPN3','CAV3','CCDC78','CFL2','CHAT','CHKB','CHRNA1','CHRNB1',
			'CHRND','CHRNE','CHRNG','CLCN1','CLN3','CNBP','CNTN1','COL6A1','COL6A2','COL6A3','COLQ','CRYAB','DAG1','DES','DMD','DMPK','DNAJB6','DNM2','DOK7','DPAGT1','DPM1','DPM2','DPM3',
			'DUX4','DYSF','EMD','FHL1','FKRP','FKTN','FLNC','GAA','GFPT1','GMPPB','GNE','HNRNPDL','HSPG2','ISCU','ISPD','ITGA7','KBTBD13','KCNA1','KCNE3','KCNJ12','KLHL40','KLHL41','KLHL9',
			'LAMA2','LAMB2','LAMP2','LARGE','LDB3','LIMS2','LMNA','LMOD3','MATR3','MCOLN1','MEGF10','MLTK','MSTN','MTM1','MTMR14','MUSK','MYH2','MYH7','MYO18B','MYOT','MYPN','NEB','ORAI1',
			'PABPN1','PLEC','POMGNT1','POMGNT2','POMK','POMT1','POMT2','PREPL','PTPLA','PTRF','RAPSN','RYR1','SCN4A','SEPN1','SGCA','SGCB','SGCD','SGCG','SMCHD1','SPEG','STAC3','STIM1','SYNE1',
			'SYNE2','SYT2','TCAP','TIA1','TMEM43','TMEM5','TNNT1','TNPO3','TOR1AIP1','TPM2','TPM3','TRAPPC11','TRIM32','TRIM54','TRIM63','TTN','VCP','VMA21']

	    f_panel_report.write('"'+'","'.join(header)+'"')
	    f_panel_report.write('\n')
	    for row in reader:
		if (row['Gene'] in muscle_genes):
		    row['UCSC_Link'] = 'UCSC_Link' #a bug with quotes
		    values = []
		    for column in header:
			values.append(row[column])
		    f_panel_report.write('"'+'","'.join(values)+'"')
		    f_panel_report.write('\n')

f_original_report.close()
f_panel_report.close()
