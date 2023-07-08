using Parameters
@with_kw mutable struct Converter

    obj_id = 8
    Conv_ID = []
    DC_Bus_ID = [] # this would be also an AC bus in case of back2back converter
    AC_Bus_ID = []
    rate = []
    GeneralSwitch = []
    type = :ACDC # can be :B2B for back to back Converter in a substation
    loss_a = 0
    loss_b = 0
    gen_dc_id = []
    gen_ac_id = []

end