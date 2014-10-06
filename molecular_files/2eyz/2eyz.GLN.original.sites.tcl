mol new 2eyz.GLN.original.xyz
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
color Name I blue2
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
topo addbond 0 17
topo addbond 1 0
topo addbond 2 1
topo addbond 2 18
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 6
topo addbond 9 0
topo addbond 10 1
topo addbond 11 4
topo addbond 12 4
topo addbond 13 5
topo addbond 14 5
topo addbond 15 8
topo addbond 16 8
draw arrow {444.002000 17.447000 2.662000} {444.644564 16.915472 2.110104} red
draw arrow {444.002000 17.447000 2.662000} {444.759477 17.996217 3.014972} red
draw arrow {444.002000 17.447000 2.662000} {444.117496 16.802144 3.417528} red
draw arrow {440.540000 17.575000 2.300000} {439.684606 17.707526 1.799263} blue
draw arrow {440.540000 17.575000 2.300000} {440.069887 17.782242 3.157930} blue
draw arrow {440.540000 17.575000 2.300000} {440.757472 18.544272 2.185028} blue
