sub defineGlobalVariable(oProject, varName, varValue)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:NewProps", Array("NAME:$"+varName, "PropType:=", "VariableProp", "UserDef:=",  _
  true, "Value:=", varValue))))
end sub

sub makeGBHidden(oProject, varName)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:ChangedProps", Array("NAME:$"+varName, "Hidden:=", true))))
end sub

sub uniteEdges(oEditor, eNames)

oEditor.Unite Array("NAME:Selections", "Selections:=",  _
  eNames), Array("NAME:UniteParameters", "KeepOriginals:=",  _
  false)

end sub 

sub coverLoop(oEditor, lName)
oEditor.CoverLines Array("NAME:Selections", "Selections:=", lName, "NewPartsModelFlag:=",  _
  "Model")
end sub

sub rename(oEditor, oldName, newName)
oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _
  oldName), Array("NAME:ChangedProps", Array("NAME:Name", "Value:=", newName))))
end sub

sub subtract(oEditor, toolParts, blankParts)
oEditor.Subtract Array("NAME:Selections", "Blank Parts:=", blankParts, "Tool Parts:=",  _
  toolParts), Array("NAME:SubtractParameters", "KeepOriginals:=", false)
end sub

sub changeObjectColor(oEditor, oName, R, G, B)
oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _
  oName), Array("NAME:ChangedProps", Array("NAME:Color", "R:=", R, "G:=", G, "B:=", B))))
end sub

sub drawSegment(oEditor, index1, index2, name)

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index1-1)+"]", "Y:=", "$y_pts["+cstr(index1-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "Line", "StartIndex:=", 0, "NoOfPoints:=", 2)), Array("NAME:PolylineXSection", "XSectionType:=",  _
  "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=",  _
  "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=",  _
  "Corner")), Array("NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=",  _
  "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=", "Global", "UDMId:=",  _
  "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

sub drawArcCPA(oEditor, index1, index2, index3, name)

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "1mm", "Y:=",  _
  "1mm", "Z:=", "0mm"), Array("NAME:PLPoint", "X:=",  _
  "-0.368087051817776mm", "Y:=", "-0.206384113560369mm", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "AngularArc", "StartIndex:=", 0, "NoOfPoints:=", 3, "NoOfSegments:=", "0", "ArcAngle:=",  _
  "$e_angles["+cstr(index3-1)+"]", "ArcCenterX:=", "$x_pts["+cstr(index1-1)+"]", "ArcCenterY:=", "$y_pts["+cstr(index1-1)+"]", "ArcCenterZ:=",  _
  "0mm", "ArcPlane:=", "XY")), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _
  "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _
  "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _
  name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=",  _
  "Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

sub drawArc(oEditor, index1, index2, index3, name)

arg1 = "atan2(" + "$y_pts["+cstr(index3-1)+"]" + "-" + "$y_pts["+cstr(index1-1)+"]," + "$x_pts["+cstr(index3-1)+"]" + "-" + "$x_pts["+cstr(index1-1)+"])"
arg2 = "atan2(" + "$y_pts["+cstr(index2-1)+"]" + "-" + "$y_pts["+cstr(index1-1)+"]," + "$x_pts["+cstr(index2-1)+"]" + "-" + "$x_pts["+cstr(index1-1)+"])"

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "1mm", "Y:=",  _
  "1mm", "Z:=", "0mm"), Array("NAME:PLPoint", "X:=",  _
  "-0.368087051817776mm", "Y:=", "-0.206384113560369mm", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "AngularArc", "StartIndex:=", 0, "NoOfPoints:=", 3, "NoOfSegments:=", "0", "ArcAngle:=",  _
  arg1+"-"+arg2, "ArcCenterX:=", "$x_pts["+cstr(index1-1)+"]", "ArcCenterY:=", "$y_pts["+cstr(index1-1)+"]", "ArcCenterZ:=",  _
  "0mm", "ArcPlane:=", "XY")), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _
  "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _
  "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _
  name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=",  _
  "Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

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
Set oProject = oDesktop.NewProject
oProject.InsertDesign "Maxwell 2D", "NewDesign", "Magnetostatic", ""
Set oDesign = oProject.SetActiveDesign("NewDesign")
Set oEditor = oDesign.SetActiveEditor("3D Modeler")

call defineGlobalVariable(oProject, "x_pts", "[0,28.2900600918,29.089779105,29.6421998266,42.6336904312,43.2,48,46.3644396619,27.335700884,28.3,29.8470569183,19.3,22.67852594,25.67852594,27.5,25.4066871441,17.8308749775,25.67852594,22.67852594,25.1318115072] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0,0.75,0.771201413428,3.49095919837,6.97197534562,0,0,12.4233141649,7.3245789764,0,0,0,0,0,0,10.52379439,7.38579024465,8.34345884812,8.34345884812,9.86874107313] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,-9.28748835642,0,15,0,-13.481382813,0,0,1.51861718696,6.71677945498,0,0,0,22.5,0,-22.5,0,0,0,3.43890639238,0] deg")
call makeGBHidden(oProject, "e_angles")
call drawSegment(oEditor, 2, 3, "stator_loop_1_e1")
call drawSegment(oEditor, 3, 4, "stator_loop_1_e2")
call drawSegment(oEditor, 4, 5, "stator_loop_1_e3")
call drawArcCPA(oEditor, 1, 5, 4, "stator_loop_1_e4")
call drawSegment(oEditor, 6, 7, "stator_loop_1_e5")
call drawArcCPA(oEditor, 1, 7, 6, "stator_loop_1_e6")
call drawSegment(oEditor, 8, 9, "stator_loop_1_e7")
call drawArcCPA(oEditor, 1, 9, 8, "stator_loop_1_e8")
call uniteEdges(oEditor, "stator_loop_1_e1,stator_loop_1_e2,stator_loop_1_e3,stator_loop_1_e4,stator_loop_1_e5,stator_loop_1_e6,stator_loop_1_e7,stator_loop_1_e8")
call coverLoop(oEditor, "stator_loop_1_e1")
call rename(oEditor, "stator_loop_1_e1", "stator")
call changeObjectColor(oEditor, "stator", 200,200,200)
call drawSegment(oEditor, 11, 6, "sca_loop_1_e10")
call drawArcCPA(oEditor, 1, 5, 4, "sca_loop_1_e4")
call drawSegment(oEditor, 4, 5, "sca_loop_1_e3")
call drawArcCPA(oEditor, 1, 11, 12, "sca_loop_1_e12")
call uniteEdges(oEditor, "sca_loop_1_e10,sca_loop_1_e4,sca_loop_1_e3,sca_loop_1_e12")
call coverLoop(oEditor, "sca_loop_1_e10")
call rename(oEditor, "sca_loop_1_e10", "sca")
call changeObjectColor(oEditor, "sca", 255,137,39)
call drawSegment(oEditor, 10, 11, "sap_loop_1_e9")
call drawArcCPA(oEditor, 1, 11, 12, "sap_loop_1_e12")
call drawSegment(oEditor, 3, 4, "sap_loop_1_e2")
call drawSegment(oEditor, 2, 3, "sap_loop_1_e1")
call drawArcCPA(oEditor, 1, 10, 11, "sap_loop_1_e11")
call uniteEdges(oEditor, "sap_loop_1_e9,sap_loop_1_e12,sap_loop_1_e2,sap_loop_1_e1,sap_loop_1_e11")
call coverLoop(oEditor, "sap_loop_1_e9")
call rename(oEditor, "sap_loop_1_e9", "sap")
call changeObjectColor(oEditor, "sap", 0,255,255)
call drawSegment(oEditor, 13, 14, "magnet_loop_1_e14")
call drawSegment(oEditor, 14, 18, "magnet_loop_1_e19")
call drawSegment(oEditor, 18, 19, "magnet_loop_1_e20")
call drawSegment(oEditor, 19, 13, "magnet_loop_1_e21")
call uniteEdges(oEditor, "magnet_loop_1_e14,magnet_loop_1_e19,magnet_loop_1_e20,magnet_loop_1_e21")
call coverLoop(oEditor, "magnet_loop_1_e14")
call rename(oEditor, "magnet_loop_1_e14", "magnet")
call changeObjectColor(oEditor, "magnet", 28,255,28)
call drawArcCPA(oEditor, 1, 18, 22, "rap_loop_1_e22")
call drawSegment(oEditor, 20, 19, "rap_loop_1_e23")
call drawSegment(oEditor, 18, 19, "rap_loop_1_e20")
call uniteEdges(oEditor, "rap_loop_1_e22,rap_loop_1_e23,rap_loop_1_e20")
call coverLoop(oEditor, "rap_loop_1_e22")
call rename(oEditor, "rap_loop_1_e22", "rap")
call changeObjectColor(oEditor, "rap", 0,255,255)
call drawSegment(oEditor, 12, 13, "rotor_loop_1_e13")
call drawSegment(oEditor, 19, 13, "rotor_loop_1_e21")
call drawSegment(oEditor, 20, 19, "rotor_loop_1_e23")
call drawArcCPA(oEditor, 1, 18, 22, "rotor_loop_1_e22")
call drawSegment(oEditor, 14, 18, "rotor_loop_1_e19")
call drawSegment(oEditor, 14, 15, "rotor_loop_1_e15")
call drawArcCPA(oEditor, 1, 15, 16, "rotor_loop_1_e16")
call drawSegment(oEditor, 16, 17, "rotor_loop_1_e17")
call drawArcCPA(oEditor, 1, 17, 18, "rotor_loop_1_e18")
call uniteEdges(oEditor, "rotor_loop_1_e13,rotor_loop_1_e21,rotor_loop_1_e23,rotor_loop_1_e22,rotor_loop_1_e19,rotor_loop_1_e15,rotor_loop_1_e16,rotor_loop_1_e17,rotor_loop_1_e18")
call coverLoop(oEditor, "rotor_loop_1_e13")
call rename(oEditor, "rotor_loop_1_e13", "rotor")
call changeObjectColor(oEditor, "rotor", 200,200,200)

