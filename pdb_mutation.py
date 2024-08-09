from pymol import cmd
cmd.wizard("mutagenesis")

# Set ATV sequence
cmd.fetch("C:/Users/rasmu/Desktop/summer_project/drugs/ATV/2fxe_aligned.pdb")
cmd.get_wizard().set_mode("LEU")
cmd.get_wizard().do_select("chain A and resid 10")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("LEU")
cmd.get_wizard().do_select("chain A and resid 63")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ARG")
cmd.get_wizard().do_select("chain A and resid 41")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("LEU")
cmd.get_wizard().do_select("chain A and resid 33")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ILE")
cmd.get_wizard().do_select("chain A and resid 13")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("LEU")
cmd.get_wizard().do_select("chain B and resid 10")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("LEU")
cmd.get_wizard().do_select("chain B and resid 63")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ARG")
cmd.get_wizard().do_select("chain B and resid 41")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("LEU")
cmd.get_wizard().do_select("chain B and resid 33")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ILE")
cmd.get_wizard().do_select("chain B and resid 13")
cmd.get_wizard().apply()


cmd.save("2fxe_aligned_fixed.pdb")


# Set IDV
from pymol import cmd
cmd.wizard("mutagenesis")

cmd.get_wizard().set_mode("ASN")
cmd.get_wizard().do_select("chain A and resid 37")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ASN")
cmd.get_wizard().do_select("chain B and resid 37")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("chain A and resid 7")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("chain B and resid 7")
cmd.get_wizard().apply()

cmd.save("1hsg_aligned.pdb")


# Set LPV
from pymol import cmd
cmd.wizard("mutagenesis")

cmd.get_wizard().set_mode("ASN")
cmd.get_wizard().do_select("chain A and resid 37")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ASN")
cmd.get_wizard().do_select("chain B and resid 37")
cmd.get_wizard().apply()
cmd.save("1mui_aligned.pdb")

# Set NFV
from pymol import cmd
cmd.wizard("mutagenesis")

cmd.get_wizard().set_mode("ASN")
cmd.get_wizard().do_select("chain A and resid 37")
cmd.get_wizard().apply()

cmd.get_wizard().set_mode("ASN")
cmd.get_wizard().do_select("chain B and resid 37")
cmd.get_wizard().apply()
cmd.save("1ohr_aligned.pdb")