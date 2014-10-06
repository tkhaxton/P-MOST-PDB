mol new 2jqx.LYS.roundtrip.xyz
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
topo addbond 0 22
topo addbond 1 0
topo addbond 2 1
topo addbond 2 23
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 7
topo addbond 9 0
topo addbond 10 1
topo addbond 11 4
topo addbond 12 4
topo addbond 13 5
topo addbond 14 5
topo addbond 15 6
topo addbond 16 6
topo addbond 17 7
topo addbond 18 7
topo addbond 19 8
topo addbond 20 8
topo addbond 21 8
draw arrow {10.171000 26.643000 -5.892000} {9.462877 27.253903 -5.537939} red
draw arrow {10.171000 26.643000 -5.892000} {9.562007 26.368355 -6.636109} red
draw arrow {10.171000 26.643000 -5.892000} {9.813662 25.900459 -5.325482} red
draw arrow {12.450000 27.429000 -5.059000} {11.788706 27.430974 -5.809124} blue
draw arrow {12.450000 27.429000 -5.059000} {11.913552 28.126732 -4.584243} blue
draw arrow {12.450000 27.429000 -5.059000} {12.974323 28.145357 -5.519347} blue
draw arrow {12.356000 29.847000 -4.382000} {12.291534 30.478226 -3.609084} green
draw arrow {12.356000 29.847000 -4.382000} {12.119454 29.084881 -3.779321} green
draw arrow {12.356000 29.847000 -4.382000} {13.325479 29.703022 -4.183556} green
