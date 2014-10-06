mol new 2eyz.ARG.original.xyz
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
draw arrow {463.686000 16.199000 -7.592000} {464.383109 15.532684 -7.327308} red
draw arrow {463.686000 16.199000 -7.592000} {463.910199 16.752260 -6.789729} red
draw arrow {463.686000 16.199000 -7.592000} {463.004991 15.699073 -7.056930} red
draw arrow {465.347000 17.380000 -9.134000} {465.282006 18.334564 -9.424833} blue
draw arrow {465.347000 17.380000 -9.134000} {464.349838 17.306777 -9.151488} blue
draw arrow {465.347000 17.380000 -9.134000} {465.309011 17.668871 -8.177386} blue
draw arrow {465.572000 20.064000 -11.690000} {464.877062 20.733339 -11.427234} green
draw arrow {465.572000 20.064000 -11.690000} {465.716258 20.551771 -12.550970} green
draw arrow {465.572000 20.064000 -11.690000} {464.867549 19.503586 -12.125529} green
