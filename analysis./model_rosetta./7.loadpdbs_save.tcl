proc loadpdbs_save {} {

	set path "/projects/cc/audagnotto_kblr332/sync/trRosetta/targets-processed/withdeepMSA/T1034/2.models_production/c.3dmodel-cst/b.glytoala/step_0_model_generation/"

	#read file line
	set fp [open "list_pdbs_lt-603.txt" r]
	set file_data [read $fp]
	close $fp

	set datas [split $file_data "\n"]
	puts $datas
	set npdbs [llength $datas]
	set new_npdbs [expr $npdbs -1]
	set start [lindex $datas 0]
	mol new $path/$start.pdb
	set id [molinfo top get id]
	animate delete beg 0 end 0 $id
	
	#upload the models and write the name correspondent to the frame.
	#needs later for retrieve energy
	set name_out [open "namemodel-frame.txt" w]
	for {set i 0} {$i < $new_npdbs} {incr i} {
		set pdb [lindex $datas $i]
		puts $name_out "$i $pdb"
		mol addfile $path/$pdb.pdb $id
	}
	
	
	set reference [atomselect top "name CA" frame 0]
	set selection [atomselect top "name CA" frame 0]
	set all [atomselect top all frame 0]
	set nf [molinfo top get numframes]
	
	for {set j 1} {$j < $nf} {incr j} {
		$selection frame $j
		$all frame $j
		set matrix [measure fit $selection $reference]
		$all move $matrix
	}
	animate write pdb allpdbs_lt603.pdb
	close $name_out
}
loadpdbs_save
exit
