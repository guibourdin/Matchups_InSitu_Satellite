function insitu_remote_match = Matchups_InSitu_PolySatellite(data, pathOC, varargin)
% author: Guillaume Bourdin
% created: Mar 04, 2020

% Matchups_InSitu_Satellite for polymer output
%
% INPUT:
%   - data_in: In-situ data <NxM Table> must contain at least:
%         - dt <Nx1 datenum or datetime> date and time precise to the second
%         - lat spectrum <Nx1 double> latitude in decimal degrees
%         - lon spectrum <Nx1 double> longitude in decimal degrees
%         - any other variable M x <Nx1 double>
%   - pathOC: string of path to data folder
%   - optional argument:
%       - (1) variable to extract: <1xM cellstr> exact variable matching geophysical variables
%           (default = all variables)
%       - (2) area match: size of square area around in-situ measurement (in pixel)
%           (default = 5; for 5x5 pixel area)
%       - (3) time interval before and after satellite overpass (in minutes)
%           (default = 180; for 3h)
%       - (4) plot: logical
%           (default = false)
%       - (5) variable to plot: <1x2 cellstr> string of exact in-situ data
%           variable (1) and matching OC variable (2) to plot
%            (default = plot == false)
%       - (6) image name list to match: <Nx6 struct> structure listing
%           image names to match with given by: dir('*.nc')
%            (default = all *.nc files in folder)
%
% OUTPUT:
%   - data_out: <NxM table> of matchups
% example: insitu_remote_match = Matchups_InSitu_Satellite(par, pathOC, var_to_extract, 5, 180, true, {'par', 'ipar'}, selseg);

cd(pathOC);
% list NC files in folder
if nargin < 8
    ncOCfiles = dir('*.nc');
else
    ncOCfiles = varargin{6};
end
NOCfiles = length(ncOCfiles); 

if nargin < 1
    error('missing in-situ data & image list')
elseif nargin < 2
    error('missing path to image')
elseif nargin < 3
    varargin{2} = 5;
    varargin{3} = 180;
    varargin{4} = false;
    NCinf = ncinfo(ncOCfiles(1).name);
    varargin{1} = fullfile({NCinf.Variables.Name});
    varargin{1}(contains(varargin{1}, {'latitude', 'longitude'})) = [];
    warning('missing variable to extract; default = all variables')
    warning('missing matching area; default = 5x5 pixels')
    warning('missing matching time interval; default = 3h')
    warning('plot set to false')
elseif nargin < 4
    varargin{2} = 5;
    varargin{3} = 180;
    varargin{4} = false;
    warning('missing matching area; default = 5x5 pixels')
    warning('missing matching time interval; default = 180 minutes')
    warning('plot set to false')
elseif nargin < 5
    varargin{3} = 180;
    varargin{4} = false;
    warning('missing matching time interval; default = 180 minutes')
    warning('plot set to false')
elseif nargin < 6
    varargin{4} = false;
    warning('plot set to false')
elseif all(any(nargin < 7 | isempty(varargin{5})) & varargin{4})
    varargin{4} = false;
    warning('missing variable to plot: plot set to false')
elseif all(nargin > 6 & ~isempty(varargin{5}))
    varargin{5}{1} = ['insitu_' varargin{5}{1}];
    varargin{5}{2} = ['remote_' varargin{5}{2}];
end

if isempty(varargin{1})
    warning('missing variable to extract; default = all variables')
    NCinf = ncinfo(ncOCfiles(1).name);
    varargin{1} = fullfile({NCinf.Variables.Name});
    varargin{1}(contains(varargin{1}, {'latitude', 'longitude'})) = [];
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
unit = cell(1, size(data,2)+5+size(varargin{1},2));
delete(gcp('nocreate'))
local_pool = parpool;
parfor(i = 1:NOCfiles, local_pool.NumWorkers)
% for i = 1:NOCfiles
    cd(pathOC);
    fprintf('Processing %s', ncOCfiles(i).name)
    try 
        if contains(ncOCfiles(i).name, 'L2')
            NCinf = ncinfo(ncOCfiles(i).name);
            att = {NCinf.Attributes.Value}';
            if any(contains({NCinf.Attributes.Name}', {'start_time', 'stop_time'}))
                datetimeOC = mean([datetime(att{contains({NCinf.Attributes.Name}', 'start_time')}) ...
                    datetime(att{contains({NCinf.Attributes.Name}', 'stop_time')})]);
            elseif any(contains({NCinf.Attributes.Name}', {'sensing_time'}))
                datetimeOC = datetime(att{contains({NCinf.Attributes.Name}', 'sensing_time')});
            end
        elseif contains(ncOCfiles(i).name, 'L3m')
            error('L3 Polymer matchup not implemented yet')
            % exatract lat/lon from L3 images
%             datetimeOC = (datetime(ncreadatt(ncOCfiles(i).name, '/', 'time_coverage_start'), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z'):...
%                 minutes(1):...
%                 datetime(ncreadatt(ncOCfiles(i).name, '/', 'time_coverage_end'), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z'))';
        end
        % identify insitu lat lon for satellite overpass time
        t_match = find(data.dt > datetimeOC(1) - varargin{3} & data.dt < datetimeOC(end) + varargin{3});
        sel_datalon = data.lon(t_match);
        if any(t_match) % if time match
            fprintf(' - time matched on %s', datestr(median(data.dt(t_match)),'yyyy/mm/dd'))
            if contains(ncOCfiles(i).name, 'L2')
                % exatract lat/lon from L2 image
                latOCini = ncread(ncOCfiles(i).name,'/latitude');
                lonOCini = ncread(ncOCfiles(i).name,'/longitude');
            elseif contains(ncOCfiles(i).name, 'L3m')
                % exatract lat/lon from L3 images
                latv = ncread(ncOCfiles(i).name, '/lat');
                lonv = ncread(ncOCfiles(i).name, '/lon');
                lonOCini = ones(length(latv),1).*lonv';
                latOCini = (latv.*ones(1,length(lonv)));
            end
            latOC = latOCini(:);
            lonOC = lonOCini(:);
            sellonOC = lonOCini;
            sellonOC(sellonOC<0) = sellonOC(sellonOC<0) + 360;
            sel_datalon(sel_datalon<0) = sel_datalon(sel_datalon<0) + 360;
            % check if tara is in the image
            insitu_in_image = inpolygon(sel_datalon, data.lat(t_match), sellonOC(:), latOC);
            if any(insitu_in_image) % if all insitu match LatLon are in image
                fprintf(' Lat/Lon matched')
                % remove tmatch out of image
                t_match(~insitu_in_image) = [];
                raddisID = acos(sin(median(data.lat(t_match))*pi/180).*sin(latOC*pi/180)...
                    +cos(median(data.lat(t_match))*pi/180).*cos(latOC*pi/180)...
                    .*cos(abs(median(data.lon(t_match))*pi/180-lonOC*pi/180)));%calculate the distance in radians with tsg at for each iteration
                [~,n_indices] = sort(raddisID(:),'ascend');    % sort the distances to find closest match area
                K = n_indices(1:varargin{2}^2,:);
                kmdis = rad2km(raddisID(K));% * 3437.74677 (nautical miles); % calculate the closest distance in km
                reso_sat = rad2km(acos(sin(latOCini(K(1))*pi/180).*sin(latOCini(K(1)+1)*pi/180)...
                    +cos(latOCini(K(1))*pi/180).*cos(latOCini(K(1)+1)*pi/180)...
                    .*cos(abs(lonOCini(K(1))*pi/180-lonOCini(K(1)+1)*pi/180))));%calculate the distance in km
                if min(kmdis) < 1.2*reso_sat % matchup if closest match is closer than 1.2 * satellite resolution
                    temp = cell(1, size(data,2)+5+size(varargin{1},2));
                    for j = 1:size(data,2)
                        temp{1,j} = data.(data.Properties.VariableNames{j})(t_match);
                    end
                    remote = cell(size(varargin{1},2),1);
                    unit = cell(1, size(data,2)+5+size(varargin{1},2));
                    for k = 1:size(varargin{1},2)
                        if contains(ncOCfiles(i).name, 'L2')
                            % exatract data from L2 image
                            remote{k} = ncread(ncOCfiles(i).name,['/' varargin{1}{k}]);
                        elseif contains(ncOCfiles(i).name, 'L3m')
                            % exatract data from L3 images
                            remote{k} = ncread(ncOCfiles(i).name, ['/' varargin{1}{k}]);
                        end                     
                        try
                            unit{size(data,2)+4+k} = strrep(ncreadatt(ncOCfiles(i).name,['/' varargin{1}{k}], 'description'), '^-', '^-^');
                        catch
                        end
                        temp{1,size(data,2)+4+k} = remote{k}(K);
                    end
                    temp{1,size(data,2)+1} = ncOCfiles(i).name;
                    temp{1,size(data,2)+2} = datetimeOC;
                    temp{1,size(data,2)+3} = latOC(K);
                    temp{1,size(data,2)+4} = lonOC(K);
                    temp{1,end} = kmdis;

                    insitu_remote_match(i,:) = temp;
                    
                    fprintf(' done\n')

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
                            idt = [strcmp(data.Properties.VariableNames,'dt') false(1,size(temp,2)-size(data,2))];
                            ilat = [strcmp(data.Properties.VariableNames,'lat') false(1,size(temp,2)-size(data,2))];
                            ilon = [strcmp(data.Properties.VariableNames,'lon') false(1,size(temp,2)-size(data,2))];
                            temp{1,size(data,2)+4}(temp{1,size(data,2)+4}<0) = temp{1,size(data,2)+4}...
                                (temp{1,size(data,2)+4}<0)+360;
                            temp{1,ilon}(temp{1,ilon}<0)...
                                = temp{1,ilon}(temp{1,ilon}<0)+360;
                            colorbar(ax2,'northoutside','Visible','off'); cb2 = colorbar(ax2,'eastoutside','FontSize',15);
                            %%%%%%%%%%% change marker size
                            % plot area of match
    %                         m_scatter(temp{1,size(data,2)+4},temp{1,size(data,2)+3},...
    %                             100,[0.5 0.5 0.5],'filled', 'Marker', 's');
                            pat = [temp{1,size(data,2)+4}(temp{1,end}>median(temp{1,end}))...
                                temp{1,size(data,2)+3}(temp{1,end}>median(temp{1,end}))];
                            c = mean(pat,1);
                            d = pat-c ;
                            th = atan2(d(:,2),d(:,1));
                            [~, idpat] = sort(th);
                            pat = pat(idpat,:);
                            patf = [pat; pat(1,:)];
                            m_patch(patf(:,1),patf(:,2), [0.85 0.85 0.85],'FaceAlpha',0.6);
                            ylabel(cb2, strrep(sprintf('%s (%s)', varargin{5}{1}, unit{id_plotinsitu}),'_',' '));
                            if isempty(unit{id_plotinsitu})
                                warning('In-situ unit missing')
                            end
                            %%%%%%%%%%% change marker size
                            m_gshhs_i('patch',[0.7 0.7 0.7],'edgecolor',[0.7 0.7 0.7],'parent',ax);
                            sel_match_dt = temp{1,idt} > datetimeOC(1) & temp{1,idt} < datetimeOC(end);
                            m_scatter(temp{1,ilon}(sel_match_dt),temp{1,ilat}(sel_match_dt),...
                                100,temp{1,id_plotinsitu}(sel_match_dt),'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                            m_scatter(temp{1,ilon}(sel_match_dt),temp{1,ilat}(sel_match_dt),...
                                100,temp{1,id_plotinsitu}(sel_match_dt),'filled');
    %                         caxis(ax, [min(min(remote{id_plotremote(size(data,2)+5:end-1)}))...
    %                             max(max(remote{id_plotremote(size(data,2)+5:end-1)}))]);
    %                         caxis(ax2, [min(temp{1,id_plotinsitu})...
    %                             max(temp{1,id_plotinsitu})]);
                            selday_data = data.dt > median(temp{1,idt}) - days(0.5)...
                                & data.dt < median(temp{1,idt}) + days(0.5);
                            fhv=figure(2);
                            set(fhv,'units','normalized','outerposition',[0.5 0.025 0.5 0.975]);
                            scatter(data.dt(selday_data),data.(varargin{5}{1}(8:end))(selday_data),7,'filled'); hold on
                            scatter(temp{1,idt}(sel_match_dt),temp{1,id_plotinsitu}(sel_match_dt),10,'r','filled');
                            y_lim = ylim();
                            area([min(temp{1,size(data,2)+2}) max(temp{1,size(data,2)+2})],...
                                [y_lim(2), y_lim(2)], y_lim(1),'FaceColor', [0.9 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                            ylabel(strrep(sprintf('%s (%s)', varargin{5}{1}, unit{id_plotinsitu}),'_',' '));
                            if isempty(unit{id_plotinsitu})
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
                    fprintf(' no LatLon match lat=%f lon=%f\n', ...
                        median(median(latOC)), median(median(lonOC)))
                end
            else
                fprintf(' no LatLon match lat=%f lon=%f\n', ...
                        median(median(latOC)), median(median(lonOC)))
            end
        else
            fprintf(' no time match\n')
        end
    catch
        fprintf(' - corrupted\n')
    end
%     clear functions
%     clear latOCini lonOCini sellonOC latOC lonOC
    latOCini = [];
    lonOCini =[];
    sellonOC = [];
    latOC = [];
    lonOC = [];
end
if exist('local_pool', 'var')
    delete(local_pool)
end
insitu_remote_match(all(cellfun(@isempty,insitu_remote_match),2),:) = [];
insitu_remote_match = cell2table(insitu_remote_match,'VariableNames', varnames);
if exist('unit','var')
    unit(cellfun(@isempty,unit))={'na'};
    insitu_remote_match.Properties.VariableUnits = unit;
end