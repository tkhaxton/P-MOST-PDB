mol new 2eyz.GLU.roundtrip.xyz
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
color Name I blue2
color Name L blue3
color Name P blue
color Name X black
proc vmd_draw_arrow {mol start end color} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.800000 [vecsub $end $start]]]
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius 0.100000
    graphics $mol cone $middle $end radius 0.200000
}
topo addbond 0 15
topo addbond 1 0
topo addbond 2 1
topo addbond 2 16
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 5
topo addbond 7 6
topo addbond 8 6
topo addbond 9 0
topo addbond 10 1
topo addbond 11 4
topo addbond 12 4
topo addbond 13 5
topo addbond 14 5
draw arrow {466.482000 18.522000 -3.518000} {465.693822 18.460567 -4.130374} red
draw arrow {466.482000 18.522000 -3.518000} {466.513421 19.511691 -3.657727} red
draw arrow {466.482000 18.522000 -3.518000} {467.096645 18.392629 -4.296123} red
draw arrow {464.273000 18.351000 -0.495000} {464.240163 18.389443 0.503721} blue
draw arrow {464.273000 18.351000 -0.495000} {463.389372 17.882940 -0.506036} blue
draw arrow {464.273000 18.351000 -0.495000} {464.740037 17.468140 -0.445661} blue
