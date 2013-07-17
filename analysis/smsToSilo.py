# Run 'visit -cli' and then 'execfile("smsToSilo.py")'
import os
import glob
smsFiles = glob.glob('/home/chaako/transferEPS/densityGradient33/densityGradient/*.sms');
print smsFiles
for smsFile in smsFiles:
    OpenDatabase(smsFile)
    AddPlot("Pseudocolor","ionDensity")
    DrawPlots()
    e = ExportDBAttributes()
    e.db_type = "Silo"
    e.variables = ("potential", "ionDensity", "eFieldX", "eFieldY", "eFieldZ",
                   "ionVelocityX", "ionVelocityY", "ionVelocityZ",
                   "ionTemperature", "referenceElectronDensity")
    #e.variables = ("potential", "ionDensity", "electronDensity", "density")
    e.filename = os.path.splitext(smsFile)[0]
    ExportDatabase(e)
    e.db_type = "VTK"
    ExportDatabase(e)

# # Database can then be opened with
# OpenDatabase("fromLoki/paddle19/ang=0.125/paddle_iter*.silo database")
# numberOfStates = TimeSliderGetNStates()
# potRefExpression = "pos_cmfe(<[" + str(numberOfStates-1) + "]i:potential>, mesh, 0.000000)"
# DefineScalarExpression("potRef", potRefExpression)
# DefineScalarExpression("dPot", "potential-potRef")
# AddPlot("Pseudocolor","dPot")
# DrawPlots()
# dPotMax = []
# for state in range(numberOfStates):
#     SetTimeSliderState(state)
#     Query("MinMax")
#     minMax = GetQueryOutputValue()
#     dPotMax = dPotMax + [max(abs(minMax[0]),abs(minMax[1]))]

# print dPotMax

# AddPlot("Pseudocolor","ionDensity")
# DrawPlots()
# ionDensities = []
# numberOfStates = TimeSliderGetNStates()
# #numberOfStates = 1
# for state in range(numberOfStates):
#     print state
#     ionDensities = ionDensities + [[]]
#     SetTimeSliderState(state)
#     Query("NumNodes")
#     numberOfNodes = int(GetQueryOutputValue())
# #    numberOfNodes = 100
#     for node in range(numberOfNodes):
#         ionDensities[state] = ionDensities[state] + [node]
#         pick = PickByNode(node)
#         ionDensities[state][node] = pick['ionDensity']
