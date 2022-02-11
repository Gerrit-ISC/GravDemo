function cb_switch(btn)
%%cb_switch (DemonstratorSimscape_vXX.slx)


current_color = get_param(btn,'BackgroundColor');
if strcmp(current_color,'green')
	%return
end

btn_tag = get_param(btn,'Tag');
sim_model = get_param(btn,'Parent');
btn_name = get_param(btn,'Name');


switch btn_tag
	
	case 'btn_load' % Single | Double
		
		hSingle = [sim_model '/Transmission System/Anchor and Load/Load/SINGLE'];
		hDouble = [sim_model '/Transmission System/Anchor and Load/Load/DOUBLE'];
		hPulley4 = [sim_model '/Transmission System/Roof Pulley Set/Sheave_04'];
		hEnds = [sim_model '/Transmission System/Roof Pulley Set/ENDS'];
		hCable = [sim_model '/Transmission System/Roof Pulley Set/Cable Winch'];
		
		if strcmp(btn_name,'Single')
			
			set_param(hEnds,'commented','on')
			set_param(hPulley4,'commented','off')
			set_param(hCable,'commented','on')
			set_param(hSingle,'commented','off')
			set_param(hDouble,'commented','on')
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Double'],'BackgroundColor','white');
			
		else%strcmp(btn_name,'Double')
			
			set_param(hEnds,'commented','off')
			set_param(hPulley4,'commented','on')
			set_param(hCable,'commented','off')
			set_param(hSingle,'commented','on')
			set_param(hDouble,'commented','off')
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Single'],'BackgroundColor','white');
			
		end
		
	case 'btn_cableflex' % Rigid | Flexi1
		
		sys = [sim_model '/Transmission System/Anchor and Load'];
		
		hDisp_L = [sim_model '/Transmission System/Anchor and Load/FlexDisp_Left'];
		hDisp_R = [sim_model '/Transmission System/Anchor and Load/FlexDisp_Right'];
		hFlex_L = [sim_model '/Transmission System/Anchor and Load/CableFlex_Left'];
		hFlex_R = [sim_model '/Transmission System/Anchor and Load/CableFlex_Right'];
		hRigid_L = [sim_model '/Transmission System/Anchor and Load/CableRigid_Left'];
		hRigid_R = [sim_model '/Transmission System/Anchor and Load/CableRigid_Right'];
		
		if strcmp(btn_name,'Rigid')
			
			set_param(hDisp_L,'commented','on')
			set_param(hDisp_R,'commented','on')
			set_param(hFlex_L,'commented','on')
			set_param(hFlex_R,'commented','on')
			set_param(hRigid_L,'commented','off')
			set_param(hRigid_R,'commented','off')

			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Flexi1'],'BackgroundColor','white');
			
		else %strcmp(btn_name,'Flexi1')
			
			set_param(hDisp_L,'commented','off')
			set_param(hDisp_R,'commented','off')
			set_param(hFlex_L,'commented','off')
			set_param(hFlex_R,'commented','off')
			set_param(hRigid_L,'commented','on')
			set_param(hRigid_R,'commented','on')

			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Rigid'],'BackgroundColor','white');

		end
		
	case 'btn_cabledamping' % cd_Const | cd_FreqDep
		
		hSW = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','SW_b');
		
		if strcmp(btn_name,'cd_Const')
			
			for i=1:length(hSW)
				set_param(hSW{i},'sw','0')
			end
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/cd_FreqDep'],'BackgroundColor','white');
			
		else %strcmp(btn_name,'cd_FreqDep')
			
			for i=1:length(hSW)
				set_param(hSW{i},'sw','1')
			end
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/cd_Const'],'BackgroundColor','white');
			
		end
				
	case 'btn_cablemass' % cm_No | cm_Yes
		
		hCM = [sim_model '/Transmission System/Anchor and Load/cable mass'];
		hC = [sim_model '/Transmission System/Anchor and Load/Load/SINGLE/Cable'];
		hCL = [sim_model '/Transmission System/Anchor and Load/Load/DOUBLE/Cable_Left'];
		hCR = [sim_model '/Transmission System/Anchor and Load/Load/DOUBLE/Cable_Right'];
		
		if strcmp(btn_name,'cm_No')
			
			set_param(hCM,'commented','on')
			set_param(hC,'commented','on')
			set_param(hCL,'commented','on')
			set_param(hCR,'commented','on')
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/cm_Yes'],'BackgroundColor','white');
			
		else%strcmp(btn_name,'cm_Yes')
			
			set_param(hCM,'commented','off')
			set_param(hC,'commented','off')
			set_param(hCL,'commented','off')
			set_param(hCR,'commented','off')
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/cm_No'],'BackgroundColor','white');
			
		end
		
	case 'btn_winch' % w_Const | w_Var
		
		hL = [sim_model '/EM System/Winch Left (yellow)/Variable Winch Inertia'];
		hR = [sim_model '/EM System/Winch Right (blue)/Variable Winch Inertia'];
		hswL = [sim_model '/EM System/Motor Winch Left/SW_wr'];
		hswR = [sim_model '/EM System/Motor Winch Right/SW_wr'];
		
		if strcmp(btn_name,'w_Const')
			
			set_param(hL,'commented','on')
			set_param(hR,'commented','on')
			% fix winch radius to const since cannot change spool
            %set_param(hswL,'sw','0') 
			%set_param(hswR,'sw','0')

			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/w_Var'],'BackgroundColor','white');
			
		else%strcmp(btn_name,'w_Var')
			
			set_param(hL,'commented','off')
			set_param(hR,'commented','off')
			%set_param(hswL,'sw','1')
			%set_param(hswR,'sw','1')

			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/w_Const'],'BackgroundColor','white');
			
		end

	case 'btn_cablelosses' % None | Combined | Individual
		
		hL = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','Bending loss');
		hSW = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','swSheavesEff');

		if strcmp(btn_name,'None')
			
			for i=1:length(hL)
				set_param(hL{i},'Commented','on')
			end
			
			for i=1:2
				set_param(hSW{i},'sw','1')
			end
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Combined'],'BackgroundColor','white');
			set_param([sim_model '/Individual'],'BackgroundColor','white');
			
		elseif strcmp(btn_name,'Combined')
			
			for i=1:length(hL)
				set_param(hL{i},'Commented','on')
			end
			
			for i=1:2
				set_param(hSW{i},'sw','0')
			end
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/None'],'BackgroundColor','white');
			set_param([sim_model '/Individual'],'BackgroundColor','white');
			
		else%strcmp(btn_name,'Individual')
			
			for i=1:length(hL)
				set_param(hL{i},'Commented','off')
			end
			
			for i=1:2
				set_param(hSW{i},'sw','1')
			end
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/None'],'BackgroundColor','white');
			set_param([sim_model '/Combined'],'BackgroundColor','white');
			
		end

	case 'btn_drive' % Simple | Huisman
		
		hS = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','Simple Drive Model');
		hH = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','PLC and Drive Logic Rev3');
		hs = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','SWTL');
		hh = find_system(sim_model,'LookUnderMasks','all','IncludeCommented','on','Name','SWTR');
		
		if strcmp(btn_name,'Simple')
			
			set_param(hS{1},'Commented','off')
			set_param(hH{1},'Commented','on')
			set_param(hs{1},'sw','1')
			set_param(hh{1},'sw','1')
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Huisman'],'BackgroundColor','white');
			
		else%strcmp(btn_name,'Huisman')

			set_param(hS{1},'Commented','on')
			set_param(hH{1},'Commented','off')
			set_param(hs{1},'sw','0')
			set_param(hh{1},'sw','0')
			
			set_param(btn,'BackgroundColor','green');
			set_param([sim_model '/Simple'],'BackgroundColor','white');
			
		end
		
		
end




