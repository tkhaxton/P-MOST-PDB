mol new 2jqx.TYR.original.xyz
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
topo addbond 0 21
topo addbond 1 0
topo addbond 2 1
topo addbond 2 22
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 5
topo addbond 8 6
topo addbond 9 7
topo addbond 10 8
topo addbond 10 9
topo addbond 11 10
topo addbond 12 0
topo addbond 13 1
topo addbond 14 4
topo addbond 15 4
topo addbond 16 6
topo addbond 17 7
topo addbond 18 8
topo addbond 19 9
topo addbond 20 11
draw arrow {26.685000 -18.942000 -2.606000} {26.392544 -18.369604 -1.839950} red
draw arrow {26.685000 -18.942000 -2.606000} {27.138136 -19.564478 -1.967887} red
draw arrow {26.685000 -18.942000 -2.606000} {27.527103 -18.408256 -2.683326} red
draw arrow {25.741000 -19.383000 -4.924000} {26.617132 -19.687992 -5.297327} blue
draw arrow {25.741000 -19.383000 -4.924000} {25.417789 -19.180060 -5.848311} blue
draw arrow {25.741000 -19.383000 -4.924000} {26.098670 -18.452519 -4.844774} blue
