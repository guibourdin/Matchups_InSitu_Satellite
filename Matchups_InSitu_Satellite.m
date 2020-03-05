function insitu_remote_match = Matchups_InSitu_Satellite(data, pathOC, varargin)
% author: Guillaume Bourdin
% created: Mar 04, 2020

% Matchups_InSitu_Satellite
%
% INPUT:
%   - data_in: In-situ data <NxM Table> must contain at least:
%         - dt <Nx1 datenum or datetime> date and time precise to the second
%         - lat spectrum <Nx1 double> latitude in decimal degrees
%         - lon spectrum <Nx1 double> longitude in decimal degrees
%         - any other variable M x <Nx1 double>
%   - pathOC: string of path to data folder
%   - optional argument: structure with three fields:
%       - (1) variable to extract: <1xM cellstr> exact variable matching OC file variables
%           (default = all OC variables)
%       - (2) area match: size of square area around in-situ measurement (in pixel)
%           (default = 5; for 5x5 pixel area)
%       - (3) time interval before and after satellite overpass (in minutes)
%           (default = 180; for 3h)
%       - (4) plot: logical
%           (default = false)
%       - (5) variable to plot: <1x2 cellstr> string of exact variable to plot
%       matching OC variable (1) and in-situ data variable (2)
%
% OUTPUT:
%   - data_out: <NxM table> of matchups
% example: insitu_remote_match = Matchups_InSitu_Satellite(par, pathOC, var_to_extract, 5, 180, true, {'par', 'ipar'});

cd(pathOC);
% list NC files in folder
ncOCfiles = dir('*.nc');
NOCfiles = length(ncOCfiles); 
NCinf = ncinfo(ncOCfiles(1).name,'/geophysical_data/');
remote_varnames = fullfile({NCinf.Variables.Name});

if nargin < 1
    error('missing in-situ data & path to OC files')
elseif nargin < 2
    error('missing path to ')
elseif nargin < 3
    varargin{2} = 5;
    varargin{3} = 180;
    varargin{4} = false;
    varargin{1} = remote_varnames;
    warning('missing variable to extract; default = all OC variables')
    warning('missing matching area; default = 5x5 pixels')
    warning('missing matching time interval; default = 3h')
    warning('plot = false')
elseif nargin < 4
    varargin{2} = 5;
    varargin{3} = 180;
    varargin{4} = false;
    warning('missing matching area; default = 5x5 pixels')
    warning('missing matching time interval; default = 180 minutes')
    warning('plot = false')
elseif nargin < 5
    varargin{3} = 180;
    varargin{4} = false;
    warning('missing matching time interval; default = 180 minutes')
    warning('plot = false')
elseif nargin < 6
    varargin{4} = false;
    warning('plot = false')
elseif nargin < 7
    warning('missing variable to plot: plot = false')
else
    varargin{5}{1} = ['insitu_' varargin{5}{1}];
    varargin{5}{2} = ['remote_' varargin{5}{2}];
end

varnames = [cellfun(@(c)['insitu_' c],data.Properties.VariableNames,'uni',false)...
    {'image_name'} {'remote_dt'} {'remote_lat'} {'remote_lon'}...
    cellfun(@(c)['remote_' c],varargin{1},'uni',false) {'dist_insitu_remote'}];

varargin{3} = minutes(varargin{3});

if ~isdatetime(data.dt)
    try
        data.dt = datetime(data.dt, 'ConvertFrom', 'datenum');
    catch
        error('date/time format not recognized')
    end
end

insitu_remote_match = cell(NOCfiles, size(data,2)+5+size(varargin{1},2));
for i = 1:NOCfiles
    cd(pathOC);
    ye = ncread(ncOCfiles(i).name,'/scan_line_attributes/year');
    D = datetime(ye,01,01)+(ncread(ncOCfiles(i).name,'/scan_line_attributes/day')-1);
    h = floor(ncread(ncOCfiles(i).name,'/scan_line_attributes/msec')/3600000);
    minu = floor(ncread(ncOCfiles(i).name,'/scan_line_attributes/msec')/60000-h*60);
    sec = floor(ncread(ncOCfiles(i).name,'/scan_line_attributes/msec')/1000)-(h*3600+minu*60);
    datetimeOC = datetime(ye,month(D),day(D),h,minu,sec);
    % identify Tara lat lon for satellite overpass time
    t_match = find(data.dt > datetimeOC(1) - varargin{3} & data.dt < datetimeOC(end) + varargin{3});
    sel_datalon = data.lon(t_match);
    if any(t_match) % if time match
        latOCini = ncread(ncOCfiles(i).name,'/navigation_data/latitude');
        lonOCini = ncread(ncOCfiles(i).name,'/navigation_data/longitude');

        latOC = latOCini(:);
        lonOC = lonOCini(:);
        sellonOC = lonOCini;
        sellonOC(sellonOC<0) = sellonOC(sellonOC<0) + 360;
        sel_datalon(sel_datalon<0) = sel_datalon(sel_datalon<0) + 360;
        % check if tara is in the image
%         ncdisp(ncOCfiles(i).name)
%         N_lat = ncreadatt(ncOCfiles(i).name,'/','northernmost_latitude');
%         S_lat = ncreadatt(ncOCfiles(i).name,'/','southernmost_latitude');
%         E_lon = ncreadatt(ncOCfiles(i).name,'/','easternmost_longitude');
%         W_lon = ncreadatt(ncOCfiles(i).name,'/','westernmost_longitude');
%         N_lat = latOC==max(latOC);
%         E_lon = sellonOC(:)==max(sellonOC(:));
%         S_lat = latOC==min(latOC);
%         W_lon = sellonOC(:)==min(sellonOC(:));
%         im_lim = [latOC(N_lat), sellonOC(N_lat); latOC(E_lon), sellonOC(E_lon);...
%             latOC(S_lat), sellonOC(S_lat); latOC(W_lon), sellonOC(W_lon);...
%             latOC(N_lat), sellonOC(N_lat);];
%         tara_in_image = inpolygon(sel_datalon,data.lat(t_match),im_lim(:,2),im_lim(:,1));
%         if all(median(data.lat(t_match)) >= min(latOC) & median(data.lat(t_match)) <= max(latOC)...
%                 & median(sel_datalon) >= min(sellonOC(:)) & median(sel_datalon) <= max(sellonOC(:))) % if LatLon match
        insitu_in_image = inpolygon([sel_datalon(1); sel_datalon(end)],...
            [data.lat(t_match(1)); data.lat(t_match(end))],sellonOC(:),latOC);
        if all(insitu_in_image) % if LatLon match
            raddisID = acos(sin(median(data.lat(t_match))*pi/180).*sin(latOC*pi/180)...
                +cos(median(data.lat(t_match))*pi/180).*cos(latOC*pi/180)...
                .*cos(abs(median(data.lon(t_match))*pi/180-lonOC*pi/180)));%calculate the distance in radians with tsg at for each iteration
            [~,n_indices]=sort(raddisID(:),'ascend');    % sort the distances to find closest match area
            K = n_indices(1:varargin{2}^2,:);
            kmdis = rad2km(raddisID(K));% * 3437.74677 (nautical miles); % calculate the closest distance in km
            if min(kmdis) < varargin{2}*2 % matchup if closest match is closer than twice the distance between center and furthest pixel
                for j = 1:size(data,2)
                    insitu_remote_match{i,j} = data.(data.Properties.VariableNames{j})(t_match);
                end
                remote = cell(size(varargin{1},2),1);
                unit = cell(1, size(data,2)+5+size(varargin{1},2));
                for k = 1:size(varargin{1},2)
                    remote{k} = ncread(ncOCfiles(i).name,['/geophysical_data/' varargin{1}{k}]);
                    try
                        unit{size(data,2)+4+k} = ncreadatt(ncOCfiles(i).name,['/geophysical_data/' varargin{1}{k}], 'units');
                    catch
                    end
                    insitu_remote_match{i,size(data,2)+4+k} = remote{k}(K);
                end
                insitu_remote_match{i,size(data,2)+1} = ncOCfiles(i).name;
                insitu_remote_match{i,size(data,2)+2} = datetimeOC;
                insitu_remote_match{i,size(data,2)+3} = latOC(K);
                insitu_remote_match{i,size(data,2)+4} = lonOC(K);
                insitu_remote_match{i,end} = kmdis;
                
                fprintf('%s - matched with %s\n', datestr(median(data.dt(t_match)),'yyyy/mm/dd'), ncOCfiles(i).name)
                if all(varargin{4} & size(varargin,2))
                    id_plotinsitu = strcmp(varnames, varargin{5}{1});
                    if all(~id_plotinsitu)
                        error('In-situ variable to plot not recognized')
                    end
                    id_plotremote = strcmp(varnames, varargin{5}{2});
                    if all(~id_plotremote)
                        error('remote variable to plot not recognized')
                    end
                    if sum(id_plotinsitu)==1
                        x = [min(sellonOC(:)) max(sellonOC(:))];
                        y = [min(latOCini(:)) max(latOCini(:))];
                        cd('C:\Users\Gui\Documents\MATLAB\Matlab_toolbox\Marc\m_map\');
                        f=figure(1);
                        set(f,'units','normalized','outerposition',[0 0.025 0.5 0.975]);
                        ax2 = axes('NextPlot','add');
                        m_proj('mercator','lon',x,'lat',y,ax2);
                        ax = axes('NextPlot','add');
                        m_proj('mercator','lon',x,'lat',y,ax);
                        axis(ax2,'off')
                        colorbar(ax,'eastoutside','Visible','off'); cb1 = colorbar(ax,'northoutside','FontSize',15);

%                         p = m_pcolor(sellonOC, latOCini, latOCini./latOCini-1); alpha(p,0.5); hold on
                        m_line([sellonOC(1,:)'; sellonOC(end,:)'; sellonOC(:,1); sellonOC(:,end)],...
                            [latOCini(1,:)'; latOCini(end,:)'; latOCini(:,1); latOCini(:,end)],...
                            'linewi',2,'color','k');     % Area outline ...
%                         m_plot(im_lim(:,2),im_lim(:,1),'Color', 'r');
                        p = m_pcolor(sellonOC, latOCini, remote{id_plotremote(size(data,2)+5:end-1)}); alpha(p,1);
                        ylabel(cb1, strrep(sprintf('%s (%s)', varargin{5}{2}, unit{id_plotremote}),'_',' '));
                        if isempty(unit{id_plotremote})
                            warning('remote unit missing')
                        end
                        m_grid('box','fancy','parent',ax);
                        m_proj('mercator','lon',x,'lat',y,ax2);
                        axis(ax2,'off');

                        set(findobj('tag','m_grid_color'),'facecolor','none');
                        colorbar(ax2,'northoutside','Visible','off');
                        colorbar(ax2,'eastoutside','Visible','off'); hold on
                        idt = [strcmp(data.Properties.VariableNames,'dt') false(1,size(insitu_remote_match,2)-size(data,2))];
                        ilat = [strcmp(data.Properties.VariableNames,'lat') false(1,size(insitu_remote_match,2)-size(data,2))];
                        ilon = [strcmp(data.Properties.VariableNames,'lon') false(1,size(insitu_remote_match,2)-size(data,2))];
                        insitu_remote_match{i,size(data,2)+4}(insitu_remote_match{i,size(data,2)+4}<0) = insitu_remote_match{i,size(data,2)+4}...
                            (insitu_remote_match{i,size(data,2)+4}<0)+360;
                        insitu_remote_match{i,ilon}(insitu_remote_match{i,ilon}<0)...
                            = insitu_remote_match{i,ilon}(insitu_remote_match{i,ilon}<0)+360;
                        colorbar(ax2,'northoutside','Visible','off'); cb2 = colorbar(ax2,'eastoutside','FontSize',15);
                        %%%%%%%%%%% change marker size
                        % plot area of match
%                         m_scatter(insitu_remote_match{i,size(data,2)+4},insitu_remote_match{i,size(data,2)+3},...
%                             100,[0.5 0.5 0.5],'filled', 'Marker', 's');
                        pat = [insitu_remote_match{i,size(data,2)+4}(insitu_remote_match{i,end}>median(insitu_remote_match{i,end}))...
                            insitu_remote_match{i,size(data,2)+3}(insitu_remote_match{i,end}>median(insitu_remote_match{i,end}))];
                        c = mean(pat,1);
                        d = pat-c ;
                        th = atan2(d(:,2),d(:,1));
                        [~, idpat] = sort(th);
                        pat = pat(idpat,:);
                        patf = [pat; pat(1,:)];
                        m_patch(patf(:,1),patf(:,2), [0.85 0.85 0.85],'FaceAlpha',0.6);
                        ylabel(cb2, strrep(sprintf('%s (%s)', varargin{5}{1}, data.Properties.VariableUnits{id_plotinsitu}),'_',' '));
                        if isempty(data.Properties.VariableUnits{id_plotinsitu})
                            warning('In-situ unit missing')
                        end
                        %%%%%%%%%%% change marker size
                        m_gshhs_i('patch',[0.7 0.7 0.7],'edgecolor',[0.7 0.7 0.7],'parent',ax);
                        sel_match_dt = insitu_remote_match{i,idt} > datetimeOC(1) & insitu_remote_match{i,idt} < datetimeOC(end);
                        m_scatter(insitu_remote_match{i,ilon}(sel_match_dt),insitu_remote_match{i,ilat}(sel_match_dt),...
                            100,insitu_remote_match{i,id_plotinsitu}(sel_match_dt),'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                        m_scatter(insitu_remote_match{i,ilon}(sel_match_dt),insitu_remote_match{i,ilat}(sel_match_dt),...
                            100,insitu_remote_match{i,id_plotinsitu}(sel_match_dt),'filled');
                        caxis(ax, [min(min(remote{id_plotremote(size(data,2)+5:end-1)}))...
                            max(max(remote{id_plotremote(size(data,2)+5:end-1)}))]);
                        caxis(ax2, [min(insitu_remote_match{i,id_plotinsitu})...
                            max(insitu_remote_match{i,id_plotinsitu})]);
                        selday_data = data.dt > median(insitu_remote_match{i,idt}) - days(0.5)...
                            & data.dt < median(insitu_remote_match{i,idt}) + days(0.5);
                        fhv=figure(2);
                        set(fhv,'units','normalized','outerposition',[0.5 0.025 0.5 0.975]);
                        scatter(data.dt(selday_data),data.par(selday_data),7,'filled'); hold on
                        scatter(insitu_remote_match{i,idt}(sel_match_dt),insitu_remote_match{i,id_plotinsitu}(sel_match_dt),10,'r','filled');
                        y_lim = ylim();
                        area([min(insitu_remote_match{i,size(data,2)+2}) max(insitu_remote_match{i,size(data,2)+2})],...
                            [y_lim(2), y_lim(2)], y_lim(1),'FaceColor', [0.9 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                        ylabel(strrep(sprintf('%s (%s)', varargin{5}{1}, data.Properties.VariableUnits{id_plotinsitu}),'_',' '));
                        if isempty(data.Properties.VariableUnits{id_plotinsitu})
                            warning('In-situ unit missing')
                        end
                        title('Validate match-up (press q) / skip matchup (press s) ');
                        [~, ~, cl] = guiSelectOnTimeSeries(fhv);
                        if any(~isempty(cl))
                            insitu_remote_match (i,:) = cell(1, size(data,2)+5+size(varargin{1},2));
                        end
                        close figure 1 figure 2
                    else
                        error('Variable to plot not recognized')
                    end
                end
            else
                fprintf('%s - No LatLon match with %s lat=%f lon=%f\n', datestr(median(data.dt(t_match)),'yyyy/mm/dd'), ...
                    ncOCfiles(i).name, median(median(latOC)), median(median(lonOC)))
            end
        else
            fprintf('%s - No LatLon match with %s lat=%f lon=%f\n', datestr(median(data.dt(t_match)),'yyyy/mm/dd'), ...
                ncOCfiles(i).name, median(median(latOC)), median(median(lonOC)))
        end
    else
        fprintf('%s - no time match\n',  datestr(median(datetimeOC),'yyyy/mm/dd'))
    end
end
insitu_remote_match(all(cellfun(@isempty,insitu_remote_match),2),:) = [];
tmp = matlab.desktop.editor.getActive; cd(fileparts(tmp.Filename));
save('insitu_remote_match.mat','insitu_remote_match')
end

