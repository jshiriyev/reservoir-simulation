from borepy.gmodel import Formation

zc = Formation()

zc["BV"] = 5600
zc["BVI"] = 5700
zc["BVII"] = 5800

zc["BV"].color = "red"
zc["BV"].hatch = ".."

zc.view()