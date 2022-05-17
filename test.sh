#Modify the first line to select the relevant files
python trajectory_converter.py -top #Topology file \
 -traj #Coordinates files \
 -r #Receptor residues -l Ligand residues \
 -ref #Reference structure
 -a #Alignment residues -ff charmm
python ifp.py
python compute_interactions.py
python training.py -f interactions_map.csv
python training.py -f interactions_map_newhyd.csv -mk mad_kernel_newhyd.sav -mm mad_ocsvm_newhyd.sav -qk qms2_kernel_newhyd.sav -qm qms2_ocsvm_newhyd.sav -r training_report_newhyd.txt
#Modify the scoring input line with the position of the file and the folder storing the docking information
python scoring.py -m *ocsvm.sav -k *kernel.sav -f /projects/cxcr4/cxcr4/classifier/dataset/adrb2/docking/functional_dataset/xqc/1048/interactions_map_m.csv -fo /projects/cxcr4/cxcr4/classifier/dataset/adrb2/docking/functional_dataset/xqc/1048/
mv MD_rescoring_0.csv MD_rescoring_0_hyd.csv
mv MD_rescoring_1.csv MD_rescoring_1_hyd.csv
python scoring.py -m *ocsvm_newhyd.sav -k *kernel_newhyd.sav -f /projects/cxcr4/cxcr4/classifier/dataset/adrb2/docking/functional_dataset/xqc/1048/interactions_map_m_newhyd.csv -fo /projects/cxcr4/cxcr4/classifier/dataset/adrb2/docking/functional_dataset/xqc/1048/ -r rescoring_report_newhyd.txt