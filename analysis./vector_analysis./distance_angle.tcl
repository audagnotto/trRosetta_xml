proc distance_angle {} {

	#output
	set outfile [open "cluster2_distance_angle.txt" w]
	#load the trajectory
	mol new cluster2.trr waitfor all
	mol addfile initial.pdb
	set id [molinfo top get id]
	set nf [molinfo $id get numframes]
	
	#selection
        set selA [atomselect $id "protein and resid 53 and name CB"]
        set selcenter1 [atomselect $id "protein and resid 53 and name NE2"]
        set selcenter2 [atomselect $id "protein and resid 130 and name NE2"]
        set selB [atomselect $id "protein and resid 130 and name CB"]

	for {set i 0} {$i < $nf} {incr i} {

		$selA frame $i
		$selcenter1 frame $i
		$selcenter2 frame $i
		$selB frame $i

		set comselA [measure center $selA]
		set comcenter1 [measure center $selcenter1]
		set comcenter2 [measure center $selcenter2]
		set comselB [measure center $selB]

		set vec1 [vecsub $comselA $comcenter1]
		set norm_vec1 [vecnorm $vec1]
		set vec2 [vecsub $comcenter2 $comselB]
		set norm_vec2 [vecnorm $vec2]
		set dot_MD [vecdot $norm_vec1 $norm_vec2]
	
		set distance [vecdist $vec1 $vec2]
		set angle [expr {acos($dot_MD) * (180/3.1415926)}]
		puts $outfile "$distance $angle $i"

		#graphics option arrow
		draw color orange
		#draw materials off
		draw material AOChalky		
		graphics $id cone $comselA $comcenter1
		graphics $id cone $comselB $comcenter2
	}
}
distance_angle
exit
