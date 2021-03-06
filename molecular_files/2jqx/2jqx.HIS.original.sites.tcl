mol new 2jqx.HIS.original.xyz
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
color Name L blue3
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
topo addbond 0 18
topo addbond 1 0
topo addbond 2 1
topo addbond 2 19
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 5
topo addbond 8 6
topo addbond 9 7
topo addbond 9 8
topo addbond 10 0
topo addbond 11 1
topo addbond 12 4
topo addbond 13 4
topo addbond 14 6
topo addbond 15 7
topo addbond 16 8
topo addbond 17 9
draw arrow {25.563000 17.291000 -2.314000} {26.139202 16.703628 -1.745682} red
draw arrow {25.563000 17.291000 -2.314000} {24.899153 17.360275 -1.569347} red
draw arrow {25.563000 17.291000 -2.314000} {25.086242 16.484653 -2.664009} red
draw arrow {26.966000 19.331000 -1.719000} {26.855671 20.323240 -1.661658} blue
draw arrow {26.966000 19.331000 -1.719000} {27.544515 19.348199 -0.903510} blue
draw arrow {26.966000 19.331000 -1.719000} {27.774176 19.454146 -2.294923} blue
