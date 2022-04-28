import os, sys

year = '2015'
month = '01'
field = 'CDF-S'
n1 = '150114'
n2 = '150115'
n3 = '150117'

f1 = open(f'{year}_{month}_{field}_{n1}.list', 'w+')
f2 = open(f'{year}_{month}_{field}_{n2}.list', 'w+')
f3 = open(f'{year}_{month}_{field}_{n3}.list', 'w+')

root_dir = f'/fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/{year}/{month}/{field}/g_band/single/lightcurves/files/'
with open(f'{year}_{month}_{field}_all.list', 'w+') as f:
	for i in os.listdir(root_dir):
		f.write(root_dir+i+'\n')
		if i.endswith(n1):
			f1.write(root_dir+i+'\n')
		elif i.endswith(n2):
			f2.write(root_dir+i+'\n')
		elif i.endswith(n3):
			f3.write(root_dir+i+'\n')
		
f1.close()
f2.close()
f3.close()
f.close()

