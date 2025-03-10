% Function to maximize the window via undocumented Java call.
% Reference: http://undocumentedmatlab.com/blog/minimize-maximize-figure-window
function MaximizeFigureWindow()
try
	% Turn off this warning "Warning: figure JavaFrame property will be obsoleted in a future release. For more information see the JavaFrame resource on the MathWorks web site.% "
	warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
	
	% 		MATLAB_Version = ver('MATLAB');
	% 		versionNumber = str2double(MATLAB_Version.Version);
	% 		if versionNumber < 8.4
	% 			FigurejFrame = get(handle(gcf),'JavaFrame');
	% 			FigurejFrame.setMaximized(true);
	% 		else
	% R2014b and later
	% Enlarge figure to full screen.
	jFrame = get(handle(gcf),'JavaFrame');
	% 			set(gcf,'Resize','off');
	drawnow;
	pause(0.1); % Give it time to render
	try
		jFrame.fHG1Client.setMaximized(true);  % HG1
	catch
		jFrame.fHG2Client.setMaximized(true);  % HG2
	end
	% 		end
catch ME
	errorMessage = sprintf('Error in program %s, function %s(), at line %d.\n\nError Message:\n%s', ...
		mfilename, ME.stack(1).name, ME.stack(1).line, ME.message);
	fprintf(1, '%s\n', errorMessage);
	WarnUser(errorMessage);
end
return; % from MaximizeFigureWindow()

%==========================================================================================================================
% Warn user via the command window and a popup message.
function WarnUser(warningMessage)
fprintf(1, '%s\n', warningMessage);
uiwait(warndlg(warningMessage));
return; % from WarnUser()
