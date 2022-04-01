#### This is a script to extract tanimoto coefficients of ~100k ligands. It uses Pandas library to extract coefficients of 
#### various filenames from all the ~100k ligand files. Input to the script is the name of the ligand and output is a pandas 
#### dataframe with all the ~100k tanimoto coefficients for that ligand.
#### Author: Dr. Swapnil Wagle, University of Michigan Ann Arbor, MI. Contact: swapnil.wagle92@gmail.com



#!/usr/bin/python3
import pandas as pd
import time
import sys, os

start_time = time.time()
in_prefix = sys.argv[1]

out_path = '/users/swagle/ligand_db/rocs_fastrocs_fr64_ge_65_out'
err_path = '/users/swagle/ligand_db/rocs_fastrocs_fr64_ge_65_err'

special_case = "1tb6_bio1"

flag_r = 0
flag_rcom = 0
flag_fr = 0
flag_fr64 = 0

rocs_paths = []
rocs_com_paths = []
fastrocs_paths = []
fastrocs_cons_paths = []

rocs_path1 = '/users/swagle/ligand_db/rocs_results/semi_final_files'

rocs_com_path1 = '/users/swagle/ligand_db/rocs_results/com_aligned'

rocs_com_path2 = '/users/swagle/ligand_db/com_aligned_rocs_results/curated_queries_output_com'

rocs_com_path31 = '/users/swagle/ligand_db/com_aligned_rocs_results/curated_db_output'
rocs_com_path32 = '/users/swagle/ligand_db/com_aligned_rocs_results'

fastrocs_path1 = '/users/swagle/ligand_db/fastrocs_results/curated_queries_txts/'

fastrocs_path11 = '/users/swagle/ligand_db/fastrocs_results/curated_db_all_queries_txt_files'
fastrocs_path12 = '/users/swagle/ligand_db/fastrocs_results/completed_3d_txt_files_queries'

fastrocs_path21 = '/users/swagle/ligand_db/fastrocs_results/curated_db_all_queries_txt_files'
fastrocs_path22 = '/users/swagle/ligand_db/fastrocs_results/missing_txts_of_3d_tanimoto_coefficients'
fastrocs_path23 = '/users/swagle/ligand_db/fastrocs_results/valid_queries_results'

fastrocs_cons_path1 = '/users/swagle/ligand_db/64_confs_fr_results/curated_queries_txts'

fastrocs_cons_path21 = '/users/swagle/ligand_db/64_confs_fr_results/results_txts_curated_db'
fastrocs_cons_path22 = '/users/swagle/ligand_db/64_confs_fr_results/tanimoto_files'

if os.path.isfile(os.path.join(rocs_path1, in_prefix + "_1.rpt")):
	flag_r = 1
	rocs_paths.append(os.path.join(rocs_path1, in_prefix + "_1.rpt"))
else:
	flag_r = 0

if os.path.isfile(os.path.join(rocs_com_path1, in_prefix + "_1.rpt")):
	flag_rcom = 1
	rocs_com_paths.append(os.path.join(rocs_com_path1, in_prefix + "_1.rpt"))
elif os.path.isfile(os.path.join(rocs_com_path2, in_prefix + "_1.rpt")):
	flag_rcom = 1
	rocs_com_paths.append(os.path.join(rocs_com_path2, in_prefix + "_1.rpt"))
elif (os.path.isfile(os.path.join(rocs_com_path31, in_prefix + "_1.rpt")) and os.path.isfile(os.path.join(rocs_com_path32, in_prefix + "_1.rpt"))):
	flag_rcom = 1
	rocs_com_paths.append(os.path.join(rocs_com_path31, in_prefix + "_1.rpt"))
	rocs_com_paths.append(os.path.join(rocs_com_path32, in_prefix + "_1.rpt"))
else:
	flag_rcom = 0

if (os.path.isfile(os.path.join(fastrocs_path1, in_prefix + ".txt"))):
	flag_fr = 1
	fastrocs_paths.append(os.path.join(fastrocs_path1, in_prefix + ".txt"))
elif (os.path.isfile(os.path.join(fastrocs_path11, in_prefix + ".txt")) and os.path.isfile(os.path.join(fastrocs_path12, in_prefix + ".txt"))):
	flag_fr = 1
	fastrocs_paths.append(os.path.join(fastrocs_path11, in_prefix + ".txt"))
	fastrocs_paths.append(os.path.join(fastrocs_path12, in_prefix + ".txt"))
elif (os.path.isfile(os.path.join(fastrocs_path21, in_prefix + ".txt")) and os.path.isfile(os.path.join(fastrocs_path22 , in_prefix + ".txt")) and os.path.isfile(os.path.join(fastrocs_path23, in_prefix + ".txt"))):
	flag_fr = 1
	fastrocs_paths.append(os.path.join(fastrocs_path21, in_prefix + ".txt"))
	fastrocs_paths.append(os.path.join(fastrocs_path22, in_prefix + ".txt"))
	fastrocs_paths.append(os.path.join(fastrocs_path23, in_prefix + ".txt"))
else:
	flag_fr = 0

if os.path.isfile(os.path.join(fastrocs_cons_path1, in_prefix + ".txt")):
	flag_fr64 = 1
	fastrocs_cons_paths.append(os.path.join(fastrocs_cons_path1, in_prefix + ".txt"))
elif (os.path.isfile(os.path.join(fastrocs_cons_path21, in_prefix + ".txt")) and os.path.isfile(os.path.join(fastrocs_cons_path22, in_prefix + ".txt"))):
	flag_fr64 = 1
	fastrocs_cons_paths.append(os.path.join(fastrocs_cons_path21, in_prefix + ".txt"))
	fastrocs_cons_paths.append(os.path.join(fastrocs_cons_path22, in_prefix + ".txt"))
else:
	flag_fr64 = 0


if (flag_r == 1 and flag_rcom == 1 and flag_fr == 1):
	out_file_path = os.path.join(out_path, in_prefix + ".out")
	err_file_path = os.path.join(err_path, in_prefix + ".err")
#	print(rocs_paths, rocs_com_paths, fastrocs_paths)
	df_ligs = pd.read_csv('/users/swagle/ligand_db/results_rearrangement/get_valid_db_indices.out', usecols = [0, 1], sep = '\s+', header = None)
	df_ligs.columns = ['idx', 'filename']
#	print(df_ligs.head(), df_ligs.dtypes, len(df_ligs.index))
	df_ligs.set_index('idx', inplace=True)

	df_rocs = (pd.read_csv(files, usecols=[0, 4], sep='\s+') for files in rocs_paths)
	concatenated_df_rocs = pd.concat(df_rocs, ignore_index = True)
	concatenated_df_rocs.columns = ['name', 'rocs']
	concatenated_df_rocs["filename"] = concatenated_df_rocs["name"].str.slice(stop = -2)
	concatenated_df_rocs.set_index('filename', inplace=True)
	del concatenated_df_rocs['name']
#	concatenated_df_rocs.to_csv(out_file_path, header=None, index=True, sep='\t', float_format='%.4f', mode='w')
	concatenated_df_rocs = concatenated_df_rocs[~concatenated_df_rocs.index.duplicated(keep='first')]
	concatenated_df_rocs.drop('Na', inplace = True)
#	concatenated_df_rocs['rocs'] = concatenated_df_rocs.rocs.astype(float)
#	concatenated_df_rocs_n = concatenated_df_rocs.convert_dtypes()
#	print(concatenated_df_rocs.head(), concatenated_df_rocs.dtypes)
#	concatenated_df_rocs.to_csv(out_file_path, header=None, index=True, sep='\t', float_format='%.4f', mode='w')

	if len(rocs_com_paths) == 1:
		df_rocs_com  = (pd.read_csv(files, usecols=[0, 4], sep='\s+') for files in rocs_com_paths)
		concatenated_df_rocs_com = pd.concat(df_rocs_com, ignore_index = True)
#		print(concatenated_df_rocs_com.to_string())
		concatenated_df_rocs_com.columns = ['name', 'rocs_com']
		concatenated_df_rocs_com["filename"] = concatenated_df_rocs_com["name"].str.slice(stop = -2)
		concatenated_df_rocs_com.set_index('filename', inplace=True)
		del concatenated_df_rocs_com['name']
		concatenated_df_rocs_com = concatenated_df_rocs_com[~concatenated_df_rocs_com.index.duplicated(keep='first')]
	#	print(concatenated_df_rocs_com.head())
	elif len(rocs_com_paths) == 2:
		df_rocs_com1 = pd.read_csv(rocs_com_paths[0] , usecols=[0, 4], sep='\s+')
		df_rocs_com1.columns = ['name', 'rocs_com']
		df_rocs_com1["filename"] = df_rocs_com1["name"].str.slice(stop = -2)
		df_rocs_com1.set_index('filename', inplace = True)
		del df_rocs_com1['name']
		df_rocs_com2 = pd.read_csv(rocs_com_paths[1] , usecols=[0, 4], sep='\s+')
		df_rocs_com2.columns = ['idx', 'rocs_com']
		df_rocs_com2["idx"] = df_rocs_com2["idx"].str.split(pat = "_", n = 1, expand = True)
		df_rocs_com2['idx'] = (df_rocs_com2['idx']).astype(int)
		df_rocs_com2.set_index('idx', inplace = True)
		df_rocs_com2 = df_rocs_com2[~df_rocs_com2.index.duplicated(keep='first')]
		df_rocs_com_conc = pd.concat([df_rocs_com2, df_ligs], axis = 1)
		df_rocs_com_conc.set_index('filename', inplace = True)
		frames = [df_rocs_com1, df_rocs_com_conc]
		concatenated_df_rocs_com = pd.concat(frames)
		concatenated_df_rocs_com = concatenated_df_rocs_com[~concatenated_df_rocs_com.index.duplicated(keep='first')]
#	print(concatenated_df_rocs_com.head(), concatenated_df_rocs_com.dtypes)

	df_fastrocs = (pd.read_csv(files, usecols=[0, 1], sep='\s+') for files in fastrocs_paths)
	concatenated_df_fastrocs = pd.concat(df_fastrocs, ignore_index = True)
	concatenated_df_fastrocs.columns = ['filename', 'fastrocs']
	concatenated_df_fastrocs.set_index('filename', inplace=True)
	concatenated_df_fastrocs = concatenated_df_fastrocs[~concatenated_df_fastrocs.index.duplicated(keep='first')]
#	print(concatenated_df_fastrocs.head(), concatenated_df_fastrocs.dtypes)
	
			#			db_name = "1tb6_bio1_GU4_GU6_GU0_GU5_GU8_GU9_GU8_GU9_GU8_GU9_GU8_GU5_GU1_GU6_GU2_GU3_0433J",  3c6s_bio4_RAM_NAG_RAM_RAM_GLC_RAM_NAG_RAM_RAM_GLC_RAM_0_0
#	df_ligs.set_index('filename', inplace=True)
	df_ligs['idx'] = df_ligs.index
	df_ligs.set_index('filename', inplace=True)
	merged_df_tanimoto = pd.concat([df_ligs, concatenated_df_rocs, concatenated_df_rocs_com, concatenated_df_fastrocs], axis = 1).reindex(df_ligs.index)
	convert_dict = {'idx':int }
	merged_df_tanimoto = merged_df_tanimoto.astype(convert_dict)
#	merged_df_tanimoto.to_csv(out_file_path, header=None, index=True, sep='\t',float_format='%.4f', mode='w')
#	print(merged_df_tanimoto.head(), merged_df_tanimoto.dtypes, len(merged_df_tanimoto.index))
#	merged_df_tanimoto.astype({'idx' :'int64'}).dtypes
#	merged_df_tanimoto_n = merged_df_tanimoto.convert_dtypes(convert_string= False)
#	print ("<{}>".format(merged_df_tanimoto.columns[1]))
	merged_df_tanimoto["max"] = merged_df_tanimoto[["rocs", "rocs_com", "fastrocs"]].max(axis=1)
	df_output = merged_df_tanimoto.loc[merged_df_tanimoto['max'] >= 0.65]
#	err_file = open (err_file_path, 'w')
#	err_file.write("ROCS\n")
#	df_output[df_output['rocs'].isna()].to_csv(err_file_path, header=None, index = True, sep = '\t', float_format='%.4f', mode='a')
#	err_file.close()
#	del df_output['rocs']
#	del df_output['rocs_com']
#	del df_output['fastrocs']
#	df_output.astype("idx":'int32', copy = False)
#	df_output['idx'] = df_output['idx'].map(lambda x: '%int' % x)

	df_output.to_csv(out_file_path, header=None, index=True, sep='\t', float_format='%.4f', mode='w')

#	print(merged_df_tanimoto.to_string())
#	print(df_output.to_string())
	print("--- %s seconds ---" % (time.time() - start_time))
else:
	err_file_path = os.path.join(err_path, in_prefix + ".err")
	err_file = open (err_file_path, 'w')
	err_file.write(rocs_paths)
	err_file.write(rocs_com_paths)
	err_file.write(fastrocs_paths)
	err_file.write(fastrocs_cons_paths)
	quit()

