mol new 2eyz.ASN.original.xyz
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
topo addbond 0 14
topo addbond 1 0
topo addbond 2 1
topo addbond 2 15
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 5
topo addbond 8 0
topo addbond 9 1
topo addbond 10 4
topo addbond 11 4
topo addbond 12 7
topo addbond 13 7
draw arrow {467.280000 30.116000 -9.216000} {467.613512 29.431284 -8.567977} red
draw arrow {467.280000 30.116000 -9.216000} {468.211280 30.462163 -9.329531} red
draw arrow {467.280000 30.116000 -9.216000} {467.133415 30.757355 -8.462888} red
draw arrow {467.513000 28.634000 -11.257000} {467.429403 27.638362 -11.215578} blue
draw arrow {467.513000 28.634000 -11.257000} {468.253419 28.544119 -11.923109} blue
draw arrow {467.513000 28.634000 -11.257000} {468.179927 28.608985 -10.512297} blue
