mol new 2eyz.TYR.original.xyz
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
draw arrow {460.053000 13.681000 1.023000} {459.211075 14.220594 1.023000} red
draw arrow {460.053000 13.681000 1.023000} {459.671920 13.086403 0.315022} red
draw arrow {460.053000 13.681000 1.023000} {459.670979 13.084936 1.729235} red
draw arrow {461.964000 15.091000 1.940000} {461.530504 15.756366 2.547758} blue
draw arrow {461.964000 15.091000 1.940000} {462.776152 15.087209 2.523434} blue
draw arrow {461.964000 15.091000 1.940000} {462.354501 15.837508 1.401265} blue
