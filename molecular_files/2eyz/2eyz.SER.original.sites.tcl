mol new 2eyz.SER.original.xyz
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
topo addbond 0 11
topo addbond 1 0
topo addbond 2 1
topo addbond 2 12
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 0
topo addbond 7 1
topo addbond 8 4
topo addbond 9 4
topo addbond 10 5
draw arrow {468.544000 21.683000 -3.308000} {468.331306 20.714355 -3.436406} red
draw arrow {468.544000 21.683000 -3.308000} {468.070186 21.900172 -4.161427} red
draw arrow {468.544000 21.683000 -3.308000} {469.398554 21.562322 -3.813148} red
draw arrow {468.932000 22.792000 -1.180000} {468.242913 22.334492 -1.742001} blue
draw arrow {468.932000 22.792000 -1.180000} {468.969309 23.544092 -1.838001} blue
draw arrow {468.932000 22.792000 -1.180000} {469.655718 22.317613 -1.681188} blue
