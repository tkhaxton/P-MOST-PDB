mol new 2jqx.ARG.original.xyz
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
color Name I blue2
color Name L blue3
color Name M blue
color Name P blue
color Name Q green2
color Name U green3
color Name V green
color Name W green
color Name X black
proc vmd_draw_arrow {mol start end color} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.800000 [vecsub $end $start]]]
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius 0.100000
    graphics $mol cone $middle $end radius 0.200000
}
topo addbond 0 24
topo addbond 1 0
topo addbond 2 1
topo addbond 2 25
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 7
topo addbond 9 8
topo addbond 10 8
topo addbond 11 0
topo addbond 12 1
topo addbond 13 4
topo addbond 14 4
topo addbond 15 5
topo addbond 16 5
topo addbond 17 6
topo addbond 18 6
topo addbond 19 7
topo addbond 20 9
topo addbond 21 9
topo addbond 22 10
topo addbond 23 10
draw arrow {24.526000 28.164000 -5.646000} {24.046289 27.346932 -5.965807} red
draw arrow {24.526000 28.164000 -5.646000} {23.717760 28.433627 -5.122504} red
draw arrow {24.526000 28.164000 -5.646000} {24.184496 28.673608 -6.435731} red
draw arrow {25.687000 27.553000 -3.434000} {25.178993 26.717056 -3.226328} blue
draw arrow {25.687000 27.553000 -3.434000} {25.584265 27.372425 -4.412181} blue
draw arrow {25.687000 27.553000 -3.434000} {26.542204 27.034742 -3.428148} blue
draw arrow {23.868000 26.568000 -0.902000} {24.747171 26.444660 -0.441733} green
draw arrow {23.868000 26.568000 -0.902000} {23.418432 26.673475 -0.015003} green
draw arrow {23.868000 26.568000 -0.902000} {23.710051 25.581257 -0.864719} green
