mol new 2eyz.Cterm.SER.original.xyz
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
topo addbond 0 12
topo addbond 1 0
topo addbond 2 1
topo addbond 3 2
topo addbond 4 1
topo addbond 5 4
topo addbond 6 2
topo addbond 7 0
topo addbond 8 1
topo addbond 9 4
topo addbond 10 4
topo addbond 11 5
draw arrow {467.015000 23.562000 9.709000} {467.565104 24.378960 9.535904} red
draw arrow {467.015000 23.562000 9.709000} {467.499749 23.080834 8.978591} red
draw arrow {467.015000 23.562000 9.709000} {466.334997 23.879893 9.048288} red
draw arrow {466.155000 23.653000 11.981000} {466.880056 23.276705 11.404203} blue
draw arrow {466.155000 23.653000 11.981000} {466.275146 24.546799 11.548925} blue
draw arrow {466.155000 23.653000 11.981000} {466.833129 23.896979 12.674265} blue
draw arrow {467.854000 24.808000 9.445000} {468.845314 24.744742 9.560306} black
draw arrow {467.854000 24.808000 9.445000} {467.945639 25.769094 9.184424} black
draw arrow {467.854000 24.808000 9.445000} {467.759663 25.076879 10.403543} black
