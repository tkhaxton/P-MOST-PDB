mol new 2jqx.PRO.roundtrip.xyz
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
topo addbond 0 14
topo addbond 1 0
topo addbond 2 1
topo addbond 2 15
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 6 0
topo addbond 7 1
topo addbond 8 4
topo addbond 9 4
topo addbond 10 5
topo addbond 11 5
topo addbond 12 6
topo addbond 13 6
draw arrow {6.936000 30.805000 5.825000} {6.254178 30.347383 6.395706} red
draw arrow {6.936000 30.805000 5.825000} {7.037121 29.973362 5.278966} red
draw arrow {6.936000 30.805000 5.825000} {7.660495 30.490412 6.438304} red
draw arrow {5.886000 31.340000 3.756000} {6.585408 30.655898 3.549039} blue
draw arrow {5.886000 31.340000 3.756000} {6.366402 31.575567 4.600821} blue
draw arrow {5.886000 31.340000 3.756000} {5.356810 30.649701 4.249402} blue
