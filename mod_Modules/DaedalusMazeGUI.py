'''
This file contain source code about the main GUI of DaedalusMaze, for Jupyter.
It allows the user to select a group of simulation modules and execute it.
It also displays the results of execution of each module at a text console.
'''

import SourceCode.DataManagement as DataManagement
import ipywidgets as widgets
import SourceCode.PanelDisplayer as PanelDisplayer 
import SourceCode.GroupsAction as GroupsAction
import os # used for checking if a file/folder exists

'''
Global Variables for this file:
'''
global CurrentGroupPanel # CurrentGroupPanel is the panel which is currently displayed. This contains visual representations of all the simulation modules inside the selected group
CurrentGroupPanel = widgets.HBox() 
global ExecutionPanelHeader # ExecutionPanelHeader contains the drop-down menu, the "Execute" button and a busy/done icon
ExecutionPanelHeader = widgets.HBox()
global ExecutionPanel # ExecutionPanel contains everything in the execute tab: the drop-down menu, the "Execute" button, the current group, the console
ExecutionPanel = widgets.VBox()
global ConsoleLabel
ConsoleLabel = widgets.Label('Console:')
global OutputTextArea # The console. Here text messages are displayed.
OutputTextArea = widgets.Textarea()
global AvailableGroups_dropdown # the drop-down menu which contains all available group names.
global BusyGIF_file # an animated gif showing busy action
BusyGIF_file = open('SourceCode/AppImages/busy.gif', 'rb').read()
global BlankGIF_file # an animated gif, just blank white
BlankGIF_file = open('SourceCode/AppImages/blank.gif', 'rb').read()
global StatusGIF # this image is a graphic representation of the status of execution. It resides right of the "Execution" button
StatusGIF = widgets.Image( value=BlankGIF_file, format='gif' )
global currentGroupImage
#currentGroupImage = widgets.Image( value=None, format='gif' )

# find out which groups of simulation modules exist. These will be displayed to user as options for execution.
# In order to do that parse the automatically fabricated GroupsAction.py file. It contains one function for each group.
AvailableGroupNames = list()
AvailableGroupNames.append(" ")
f = open("SourceCode/GroupsAction.py", "r")
for aLine in f:
    if aLine.startswith('def Execute_'):
        aGroupName = aLine[ 12 : aLine.index('(') ]
        AvailableGroupNames.append( aGroupName )
f.close()
AvailableGroupNames = list(set(AvailableGroupNames)) # remove duplicates
AvailableGroupNames.sort()



'''
This is called when the "Execute" button is pressed. 
It calls the proper function from the automatically fabricated GroupsAction.py file, according to the selected group.
This function calls all the simulation module functions of the group.
'''
def execute_button_clicked(b):
    global AvailableGroups_dropdown
    DisplayBusy()
    getattr( GroupsAction, 'Execute_'+str(AvailableGroups_dropdown.value) )()
    DisplayDone()




'''
This is called when the user selects oa group from the drop-down menu.
It calls the proper function from the automatically fabricated PanelDisplayer.py file, according to the selected group.
This function displays a visual representation of all simulation modules of a group.
'''
def AvailableGroups_dropdown_onChange(change):
    if change['type']=='change' and change['name']=='value' and ' ' not in change['new']:
        global CurrentGroupPanel
        CurrentGroupPanel.close()
        for child in CurrentGroupPanel.children: child.close()
        CurrentGroupPanel = getattr( PanelDisplayer, 'Construct_'+str(change['new']) )()
        # display an image explainig the selected group 
        groupImgFilename = "SourceCode/AppImages/"+(change['new'])+".png"
        if os.path.isfile( groupImgFilename ) == True:
            currentGroupImage = widgets.Image( value=open(groupImgFilename, 'rb').read(), format='png' )
            currentGroupImage.layout.margin = '0px 0px 15px 0px'
            ExecutionPanel.children = [ ExecutionPanelHeader, currentGroupImage, CurrentGroupPanel, ConsoleLabel, OutputTextArea ]
        else:
            ExecutionPanel.children = [ ExecutionPanelHeader, CurrentGroupPanel, ConsoleLabel, OutputTextArea ]



'''
Displays a text message at the console
'''
def Display( *arg ):
    msg = ""
    for item in arg:
        msg = msg + str(item) + " "
    OutputTextArea.value = OutputTextArea.value + msg + '\n'

'''
Displays a busy icon next to the execution button
'''
def DisplayBusy():
    global BusyGIF_file
    global StatusGIF
    StatusGIF.value = BusyGIF_file

'''
Displays an empty icon next to the execution button
'''
def DisplayDone():
    global BlankGIF_file
    global StatusGIF
    StatusGIF.value = BlankGIF_file


'''
This function bundles up all GUI elements and returns the parent of them, so that jupyter can display them.
'''	
def createAll():
    global CurrentGroupPanel 
    global ExecutionPanelHeader
    global ExecutionPanel
    global ConsoleLabel
    global OutputTextArea 
    global AvailableGroups_dropdown
    global BusyGIF_file 
    global BlankGIF_file


    # Create the Containers
    MainTab = widgets.Tab() # all visual elements are children of this one
    DataManagementPanel  = widgets.VBox()
    
    # Create the Execution Panel
    AvailableGroups_dropdown = widgets.Dropdown( # here are all module groups for execution
        description = 'Select Group:',
        options = AvailableGroupNames
    )
    ExecuteButton = widgets.Button(description='Execute', icon='gear', button_style='success')
    ExecuteButton.on_click(execute_button_clicked)
    ExecutionPanelHeader.layout.margin ='10px 10px 40px 10px' 
    ExecutionPanelHeader.children =[AvailableGroups_dropdown, ExecuteButton, StatusGIF]
    ConsoleLabel.layout.margin ='30px 0px 0px 0px' 
    OutputTextArea.layout.min_width='920px'
    OutputTextArea.layout.min_height='500px'
    OutputTextArea.layout.height='500px'
    
    # Create the Data Management panel
    DataManagementPanel.children = [ DataManagement.ConstructDataManagementPanel() ]
    
    # Connect all visual elements and fix the details
    AvailableGroups_dropdown.observe( AvailableGroups_dropdown_onChange ) # event listener
    ExecutionPanel.children += ( AvailableGroups_dropdown, )
    MainTab.children = [ExecutionPanel, DataManagementPanel ]
    MainTab.set_title(0, 'Execution')
    MainTab.set_title(1, 'Data Management')
    return MainTab
