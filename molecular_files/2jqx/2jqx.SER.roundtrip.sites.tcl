mol new 2jqx.SER.roundtrip.xyz
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
draw arrow {22.515000 30.985000 -4.025000} {23.125336 30.506382 -4.656202} red
draw arrow {22.515000 30.985000 -4.025000} {21.732517 30.744791 -4.599474} red
draw arrow {22.515000 30.985000 -4.025000} {22.638333 31.829526 -4.546118} red
draw arrow {22.430000 29.044000 -2.532000} {22.578212 30.032079 -2.573640} blue
draw arrow {22.430000 29.044000 -2.532000} {21.877801 29.161613 -1.706625} blue
draw arrow {22.430000 29.044000 -2.532000} {23.250433 28.944664 -1.968953} blue
