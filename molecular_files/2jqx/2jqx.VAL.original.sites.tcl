mol new 2jqx.VAL.original.xyz
topo clearbonds
set sel [atomselect top " name A"]
$sel set radius 1.200000
set sel [atomselect top " name B"]
$sel set radius 1.700000
set sel [atomselect top " name D"]
$sel set radius 1.200000
set sel [atomselect top " name E"]
$sel set radius 1.700000
set sel [atomselect top " name F"]
$sel set radius 1.200000
set sel [atomselect top " name G"]
$sel set radius 1.700000
set sel [atomselect top " name I"]
$sel set radius 1.200000
set sel [atomselect top " name J"]
$sel set radius 1.700000
set sel [atomselect top " name K"]
$sel set radius 1.200000
set sel [atomselect top " name L"]
$sel set radius 1.700000
set sel [atomselect top " name M"]
$sel set radius 1.200000
set sel [atomselect top " name P"]
$sel set radius 1.700000
set sel [atomselect top " name Q"]
$sel set radius 1.200000
set sel [atomselect top " name R"]
$sel set radius 1.700000
set sel [atomselect top " name T"]
$sel set radius 1.200000
set sel [atomselect top " name U"]
$sel set radius 1.700000
set sel [atomselect top " name V"]
$sel set radius 1.200000
set sel [atomselect top " name W"]
$sel set radius 1.700000
set sel [atomselect top " name X"]
$sel set radius 1.700000
color Name A red2
color Name B red2
color Name E red3
color Name F red
color Name G red
color Name M blue
color Name P blue
color Name X black
proc vmd_draw_arrow {mol start end color} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.800000 [vecsub $end $start]]]
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius 0.100000
    graphics $mol cone $middle $end radius 0.200000
}
topo addbond 0 16
topo addbond 1 0
topo addbond 2 1
topo addbond 2 17
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 4
topo addbond 7 0
topo addbond 8 1
topo addbond 9 4
topo addbond 10 5
topo addbond 11 5
topo addbond 12 5
topo addbond 13 6
topo addbond 14 6
topo addbond 15 6
draw arrow {6.992000 25.586000 -1.866000} {7.106455 26.383931 -1.274218} red
draw arrow {6.992000 25.586000 -1.866000} {6.708068 26.183123 -2.616218} red
draw arrow {6.992000 25.586000 -1.866000} {6.040011 25.503840 -1.571099} red
draw arrow {8.323000 24.805000 -1.895000} {8.540478 24.443632 -0.988294} blue
draw arrow {8.323000 24.805000 -1.895000} {8.366287 23.880545 -2.273825} blue
draw arrow {8.323000 24.805000 -1.895000} {9.298105 24.926635 -2.080406} blue
