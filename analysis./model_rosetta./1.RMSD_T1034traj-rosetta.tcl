proc RMSD_1AKEtraj-rosettadmpfold { } {

	#load the traj
	set pdbs [glob /projects/cc/audagnotto_kblr332/sync/trRosetta/targets-processed/withdeepMSA/T1034/2.models_production/c.3dmodel-cst/b.glytoala/step_0_model_generation/*.pdb]

	mol new [lindex $pdbs 0]
	set id_models [molinfo top get id]
	animate delete beg 0 end 0 $id_models
	set nf_models [llength $pdbs]
	for {set i 0} {$i < $nf_models} {incr i} {
		mol addfile [lindex $pdbs $i]
	}

	mol new /projects/cc/audagnotto_kblr332/sync/trRosetta/targets-processed/withdeepMSA/T1034/experimental/MD/6y75/3.production/renumbered-start.pdb
	mol addfile /projects/cc/audagnotto_kblr332/sync/trRosetta/targets-processed/withdeepMSA/T1034/experimental/MD/6y75/3.production/renumbered.dcd waitfor all
	set id_exp [molinfo top get id]
	set nf_models [molinfo $id_models get numframes]	
	
	#set the reference
	set reference [atomselect $id_exp "protein and resid 1 to 150 and name CA"]
	set align [atomselect $id_models "protein and resid 1 to 150 and name CA"]
	set tomove [atomselect $id_models "protein"]

	#set the select for the RMSD
	set ref_sel [atomselect $id_exp "protein and resid 1 to 150 and name CA"]
	set model_sel [atomselect $id_models "protein and resid 1 to 150 and name CA"]

	#measure the i,j rmsd
	set outfile [open "RMSD_pairframes.txt" w]
	for {set i 0} {$i < $nf_models} {incr i} {
		for {set j 0} {$j < $nf_models} {incr j} {
		#align
		
		$reference frame $j
		$align frame $i
		$tomove frame $i
		set matrix [measure fit $align $reference]
                $tomove move $matrix
		
		#measure RMSD
		$ref_sel frame $j
		$model_sel frame $i
		set RMSD [measure rmsd $ref_sel $model_sel]
		set name [lindex $pdbs $i]
		puts $outfile "$i $j $RMSD $name"

		}

	}

close $outfile
}

RMSD_1AKEtraj-rosettadmpfold
exit
