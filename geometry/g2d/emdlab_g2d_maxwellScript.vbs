sub updateGlobalVariable(oProject, varName, varValue)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:ChangedProps", Array("NAME:$"+varName, "Value:=",  _
  varValue))))
end sub

Dim oAnsoftApp
Dim oDesktop
Dim oProject
Dim oDesign
Dim oEditor
Dim oModule
Set oAnsoftApp = CreateObject("Ansoft.ElectronicsDesktop")
Set oDesktop = oAnsoftApp.GetAppDesktop()
oDesktop.RestoreWindow
Set oProject = oDesktop.GetActiveProject()


call updateGlobalVariable(oProject, "x_pts", "[0,26.4893846663,27.2890642034,27.7331468371,40.8859661033,41.4,48,46.3644396619,25.5970343967,26.5,27.8927253414,19.3,23.5143360503,25.7,24.8722539949,24.4421524688,18.3553907645,23.5143360503] mm")
call updateGlobalVariable(oProject, "y_pts", "[0,0.75,0.772641509434,2.97937803639,6.50367402337,0,0,12.4233141649,6.85870469522,0,0,0,0,0,6.47000627617,7.94173675544,5.96402799144,6.58880870233] mm")
call updateGlobalVariable(oProject, "e_angles", "[0,0,0,-9.03823700407,0,15,0,-13.3782047766,0,0,1.62179522342,6.13178082611,0,0,14.5811644112,3.41883558884,0,-18,0,0] deg")
