mol new 2eyz.PHE.roundtrip.xyz
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
topo addbond 0 20
topo addbond 1 0
topo addbond 2 1
topo addbond 2 21
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 5
topo addbond 8 6
topo addbond 9 7
topo addbond 10 8
topo addbond 10 9
topo addbond 11 0
topo addbond 12 1
topo addbond 13 4
topo addbond 14 4
topo addbond 15 6
topo addbond 16 7
topo addbond 17 8
topo addbond 18 9
topo addbond 19 10
draw arrow {469.447000 27.324000 -7.777000} {468.758162 26.688804 -7.427675} red
draw arrow {469.447000 27.324000 -7.777000} {468.733107 27.834694 -8.256112} red
draw arrow {469.447000 27.324000 -7.777000} {469.572931 26.744589 -8.582248} red
draw arrow {469.353000 28.592000 -5.569000} {468.878957 28.031374 -4.890044} blue
draw arrow {469.353000 28.592000 -5.569000} {468.958223 29.416580 -5.163761} blue
draw arrow {469.353000 28.592000 -5.569000} {468.565958 28.516064 -6.181209} blue
