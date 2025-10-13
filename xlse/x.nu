# xlse.nu - Nushell overlay for jrun

export def jrun [
  --onshellT              # Plot onshell T amplitude
  --Det:int@"nu-complete jrun-det"               # Plot determinant of I+VG (takes 0..3)
  --poles                 # Plot poles in complex energy plane
	--polesm # Plot poles in complex momentum plane
  --onshellT_single       # Plot single-channel onshell T
  --Detsing:int@"nu-complete jrun-det"               # Plot single-channel determinant
	--testwf	# Test wavefunction in momentum space
  ...rest                 # any remaining positional args
] {
  # build / run xmake first
  ^xmake

  # Reconstruct arguments to pass to julia:
  # include flags only if they were supplied.
  mut julia_args = []

  if $onshellT {
    $julia_args ++= ["--onshellT"]
  }

  if $poles {
    $julia_args ++= ["--poles"]
  }

  if $polesm {
    $julia_args ++= ["--polesm"]
  }

  if $onshellT_single {
    $julia_args ++= ["--onshellT_single"]
  }

  if ($Detsing != null) {
    $julia_args ++= ["--Detsing", ($Detsing | into string)]
  }

  if ($Det != null) {
    $julia_args ++= ["--Det", ($Det | into string)]
  }

  if not ($rest | is-empty) {
    $julia_args ++= ($rest | each {|it| $it | into string })
  }

	if $testwf { $julia_args ++= ["--testwf"] }

	if ($julia_args | is-empty) { return;  }
  # Run Julia with the constructed arguments
  ^julia cscript.jl ...$julia_args

  # Post-run: display expected images depending on which flag was used
  if $onshellT { ^kitten icat onshellT.png }
  if $poles { ^kitten icat pole.png }
  if $polesm { ^kitten icat polem.png }
  if $onshellT_single { ^kitten icat onshellT_single.png }
  if ($Detsing != null) { ^kitten icat detsing.png }

  if ($Det != null) { ^kitten icat det.png }
	if $testwf { ^kitten icat psi.png }
}

def "nu-complete jrun-det" [] {
  [
    { value: "3", description: "Sheet ++ (physical)" }
    { value: "2", description: "Sheet -+" }
    { value: "1", description: "Sheet +-" }
    { value: "0", description: "Sheet --" }
  ]
}

# -------- extern declaration for completion wiring --------
