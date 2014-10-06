mol new 2eyz.VAL.original.xyz
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
topo addbond 0 16
topo addbond 1 0
topo addbond 2 1
topo addbond 2 17
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 4
topo addbond 7 0
topo addbond 8 1
topo addbond 9 4
topo addbond 10 5
topo addbond 11 5
topo addbond 12 5
topo addbond 13 6
topo addbond 14 6
topo addbond 15 6
draw arrow {445.398000 12.692000 3.656000} {444.924120 12.361462 2.839800} red
draw arrow {445.398000 12.692000 3.656000} {445.975310 13.275283 3.084604} red
draw arrow {445.398000 12.692000 3.656000} {446.062944 11.950026 3.570417} red
draw arrow {444.378000 13.280000 4.650000} {443.999798 14.140508 4.308696} blue
draw arrow {444.378000 13.280000 4.650000} {443.454814 12.902129 4.720288} blue
draw arrow {444.378000 13.280000 4.650000} {444.309515 13.621671 5.587321} blue
