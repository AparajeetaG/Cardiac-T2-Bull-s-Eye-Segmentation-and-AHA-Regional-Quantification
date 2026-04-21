function results = make_bullseye_from_t2star_mat_AHA_aligned(matFile)
% AHA-aligned 16-segment bull's-eye workflow for T2* MAT files.
% Expected MAT contents:
%   T2_starmaps_all : [Nx Ny Nz Nt]
%   iField          : [Nx Ny Nz Nt Ncoils]
%
% Workflow:
% 1) Select one basal, one mid, one apical slice from MAG images
% 2) Select one cardiac phase
% 3) Draw endocardial and epicardial ROIs on each selected slice
% 4) Click the ANTERIOR RV insertion point on each selected slice
% 5) Create AHA-aligned segment labels (1:16)
% 6) Aggregate T2* values into a bull's-eye and save QA outputs
%
% NOTE:
% This version is designed to be anatomically aligned to the AHA 16-segment
% model, but the final orientation should still be verified visually on the
% saved QA figures for your dataset/viewer convention.

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

if nargin < 1 || isempty(matFile)
    [fn, pn] = uigetfile('*.mat', 'Select the MAT-file');
    if isequal(fn,0), error('No MAT file selected.'); end
    matFile = fullfile(pn, fn);
end

S = load(matFile);
if ~isfield(S,'T2_starmaps_all') || ~isfield(S,'iField')
    error('MAT file must contain T2_starmaps_all and iField.');
end

T2_starmaps_all = S.T2_starmaps_all;
iField = S.iField;

if ndims(T2_starmaps_all) < 4
    error('T2_starmaps_all must be 4D: [Nx Ny Nz Nt].');
end
if ndims(iField) < 5
    error('iField must be 5D: [Nx Ny Nz Nt Ncoils].');
end

[Nx, Ny, Nz, Nt] = size(T2_starmaps_all);
[~,~,Nz2,Nt2,~] = size(iField);
if Nz ~= Nz2 || Nt ~= Nt2
    error('Mismatch between T2_starmaps_all and iField dimensions.');
end

[saveDir, baseName, ~] = fileparts(matFile);
outDir = fullfile(saveDir, [baseName '_AHA_bullseye']);
if ~exist(outDir,'dir'), mkdir(outDir); end

% --- Build magnitude image stack for display ---
Mag = sqrt(sum(abs(iField).^2, 5));
magScale = max(Mag(:));
if ~isfinite(magScale) || magScale == 0, magScale = 1; end
Mag = Mag ./ magScale;
T2ms = double(T2_starmaps_all) * 1e3;
T2ms(T2ms < 0) = NaN;

%% -------------------- Slice selection on magnitude --------------------
figure('Name','Magnitude preview - phase 1','Color','w');
montage(Mag(:,:,:,1), 'DisplayRange', [0 0.8]);
colormap gray;
title('Magnitude preview (phase 1). Choose representative basal, mid, apical slice indices.');

prompt = {'Basal slice index:', 'Mid slice index:', 'Apical slice index:'};
def = {num2str(max(1, round(Nz*0.30))), num2str(max(1, round(Nz*0.50))), num2str(max(1, round(Nz*0.70)))};
answ = inputdlg(prompt, 'Select basal / mid / apical slices', 1, def);
if isempty(answ), error('Slice selection cancelled.'); end
selSlices = [str2double(answ{1}), str2double(answ{2}), str2double(answ{3})];
if any(~isfinite(selSlices)) || any(selSlices < 1) || any(selSlices > Nz) || any(mod(selSlices,1)~=0)
    error('Slice indices must be integers between 1 and %d.', Nz);
end
slice_labels = {'basal','mid','apical'};

%% -------------------- Phase selection on magnitude --------------------
midSlice = selSlices(2);
figure('Name','Magnitude preview - all phases at selected mid slice','Color','w');
montage(Mag(:,:,midSlice,:), 'DisplayRange', [0 0.8]);
colormap gray;
title(sprintf('Magnitude preview at selected mid slice %d. Choose one cardiac phase.', midSlice));

ansPhase = inputdlg({sprintf('Cardiac phase index (1..%d):', Nt)}, 'Select cardiac phase', 1, {num2str(max(1, round(Nt/2)))});
if isempty(ansPhase), error('Phase selection cancelled.'); end
phaseIdx = str2double(ansPhase{1});
if ~isfinite(phaseIdx) || phaseIdx < 1 || phaseIdx > Nt || mod(phaseIdx,1)~=0
    error('Phase index must be an integer between 1 and %d.', Nt);
end

% Reduce to 3-slice 3D volumes for the chosen phase
metricMaps  = T2ms(:,:,selSlices,phaseIdx);
displayMaps = Mag(:,:,selSlices,phaseIdx);
num_slice = numel(selSlices);

%% -------------------- ROI drawing + anterior RV insertion --------------------
Mask_myo = false(size(metricMaps));
roi_endo = cell(1,num_slice);
roi_epi  = cell(1,num_slice);
ref_pts  = nan(num_slice,2); % [x y]

for si = 1:num_slice
    imgMag = displayMaps(:,:,si);
    imgT2  = metricMaps(:,:,si);

    [BW, okMyo, endo, epi] = get_or_edit_myo_roi_nozoom_click2edit_safe(imgMag, imgT2, [], [], []);
    if ~okMyo || isempty(BW)
        error('ROI not completed for %s slice.', slice_labels{si});
    end

    [pt, okPt] = get_or_edit_point(imgMag, [NaN NaN]);
    if ~okPt || any(isnan(pt))
        error('Anterior RV insertion point not selected for %s slice.', slice_labels{si});
    end

    Mask_myo(:,:,si) = logical(BW);
    roi_endo{si} = endo;
    roi_epi{si}  = epi;
    ref_pts(si,:) = pt;
end

%% -------------------- Compute AHA reference angles --------------------
RVrefAngleDeg = zeros(1, num_slice);
for si = 1:num_slice
    BWmyo = Mask_myo(:,:,si);
    C = largest_component_props(BWmyo);
    cx = C.Centroid(1);
    cy = C.Centroid(2);
    x_ins = ref_pts(si,1);
    y_ins = ref_pts(si,2);
    RVrefAngleDeg(si) = mod(atan2(-(y_ins - cy), (x_ins - cx)) * 180/pi, 360);
end

%% -------------------- AHA segmentation --------------------
optsAHA = struct();
optsAHA.clockwise = false;     % anticlockwise traversal to match AHA display convention more closely
optsAHA.validRange = [5 100];  % T2* ms range
optsAHA.zero_tol = 1e-12;
optsAHA.frac_bad_threshold = 0.30;
optsAHA.abs_bad_threshold = 3;

[Segmentpix, stats_local, Mask_index_local, Mask_index_global, qa_info, report_file] = ...
    AHASegmentation_t2star_invivo_AHAcorrect(metricMaps, Mask_myo, RVrefAngleDeg, ...
    slice_labels, baseName, outDir, phaseIdx, optsAHA);

segLabelMaps = double(Mask_index_global);

%% -------------------- Original-style bull's-eye aggregation --------------------
summaryFcn = @(x) median(x,'omitnan');
out = compute_bullseye_by_ring(metricMaps, segLabelMaps, 1, 2, 3, summaryFcn);

T = table((1:16)', out.all16(:), out.base16(:), out.mid16(:), out.apex16(:), ...
    'VariableNames', {'Segment','AllSlices','Base','Mid','Apex'});

%% -------------------- QA figures --------------------
% Segment overlay figure
figQA = figure('Color','w','Name','AHA segment QA');
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','loose');
for si = 1:num_slice
    ax = nexttile(tl, si);
    imagesc(ax, displayMaps(:,:,si), [0 0.8]);
    axis(ax,'image'); axis(ax,'off'); colormap(ax,'gray'); hold(ax,'on');
    visboundaries(ax, Mask_myo(:,:,si), 'Color', 'y', 'LineWidth', 0.8);
    h = imagesc(ax, Mask_index_global(:,:,si));
    set(h,'AlphaData', 0.25*(Mask_index_global(:,:,si)>0));
    colormap(ax, gray);
    plot(ax, ref_pts(si,1), ref_pts(si,2), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
    overlay_aha_labels(ax, Mask_index_global(:,:,si), slice_labels{si});
    title(ax, sprintf('%s (orig slice %d)', upper(slice_labels{si}), selSlices(si)), 'FontSize', 11, 'FontWeight', 'bold');
end
annotation(figQA, 'textbox', [0.25 0.955 0.50 0.04], 'String', sprintf('AHA QA overlay | phase %d', phaseIdx), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 18);
annotation(figQA, 'textbox', [0.18 0.005 0.64 0.03], 'String', 'Magnitude images with myocardial ROI overlay and AHA segment labels', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11);
save_current_figure(fullfile(outDir, [baseName '_AHA_QA_overlay.png']));

% Bulls-eye figure
figBull = figure('Color','w','Name','Bull''s-eye');
plot_bullseye_aha16_AHA(out.all16, 'T2* bull''s-eye (AHA-aligned, 16 segments)');
save_current_figure(fullfile(outDir, [baseName '_bullseye_AHA16.png']));

%% -------------------- Save outputs --------------------
writetable(T, fullfile(outDir, [baseName '_bullseye_table.csv']));
save(fullfile(outDir, [baseName '_bullseye_input.mat']), 'metricMaps', 'segLabelMaps');
save(fullfile(outDir, [baseName '_bullseye_results.mat']), ...
    'T','out','metricMaps','segLabelMaps','Mask_myo','Mask_index_local','Mask_index_global', ...
    'roi_endo','roi_epi','ref_pts','RVrefAngleDeg','stats_local','Segmentpix', ...
    'qa_info','report_file','selSlices','phaseIdx','slice_labels','optsAHA');

% Console checklist
fprintf('\n===== RUN CHECKLIST =====\n');
fprintf('Selected slices: basal=%d, mid=%d, apical=%d\n', selSlices(1), selSlices(2), selSlices(3));
fprintf('Selected phase: %d\n', phaseIdx);
fprintf('Saved QA overlay: %s\n', fullfile(outDir, [baseName '_AHA_QA_overlay.png']));
fprintf('Saved bull''s-eye: %s\n', fullfile(outDir, [baseName '_bullseye_AHA16.png']));
fprintf('Saved table: %s\n', fullfile(outDir, [baseName '_bullseye_table.csv']));
fprintf('Saved input MAT for original-style aggregation: %s\n', fullfile(outDir, [baseName '_bullseye_input.mat']));

results = struct();
results.T = T;
results.out = out;
results.metricMaps = metricMaps;
results.segLabelMaps = segLabelMaps;
results.Mask_myo = Mask_myo;
results.Mask_index_local = Mask_index_local;
results.Mask_index_global = Mask_index_global;
results.stats_local = stats_local;
results.Segmentpix = Segmentpix;
results.ref_pts = ref_pts;
results.RVrefAngleDeg = RVrefAngleDeg;
results.selSlices = selSlices;
results.phaseIdx = phaseIdx;
results.slice_labels = slice_labels;
results.outDir = outDir;
results.qa_info = qa_info;
results.report_file = report_file;
end

%% ========================= LOCAL FUNCTIONS =========================

function [Segmentpix, stats, Mask_index_local, Mask_index_global, qa_info, report_file] = ...
    AHASegmentation_t2star_invivo_AHAcorrect(Imgin, Maskin, RVrefAngleDeg, slice_labels, filename, save_path, ph, opts)

    if nargin < 8 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'clockwise'), opts.clockwise = true; end
    if ~isfield(opts,'validRange'), opts.validRange = [5 100]; end
    if ~isfield(opts,'zero_tol'), opts.zero_tol = 1e-12; end
    if ~isfield(opts,'frac_bad_threshold'), opts.frac_bad_threshold = 0.30; end
    if ~isfield(opts,'abs_bad_threshold'), opts.abs_bad_threshold = 3; end

    nSlices = size(Imgin,3);
    Mask_index_local  = zeros(size(Maskin));
    Mask_index_global = zeros(size(Maskin));
    Segmentpix = cell(16, nSlices);
    stats = nan(5,16);
    report_file = '';

    qa_info = struct();
    qa_info.n_total_raw   = zeros(16, nSlices);
    qa_info.n_zero        = zeros(16, nSlices);
    qa_info.n_nan         = zeros(16, nSlices);
    qa_info.n_inf         = zeros(16, nSlices);
    qa_info.n_nonfinite   = zeros(16, nSlices);
    qa_info.n_valid_clean = zeros(16, nSlices);
    qa_info.frac_bad      = zeros(16, nSlices);
    qa_info.flag_bad      = false(16, nSlices);

    for m = 1:nSlices
        Img  = double(Imgin(:,:,m));
        Mask = logical(Maskin(:,:,m));
        if ~any(Mask(:)), continue; end

        [Y, X] = ndgrid(1:size(Mask,1), 1:size(Mask,2));
        x_c = mean(X(Mask), 'omitnan');
        y_c = mean(Y(Mask), 'omitnan');

        theta = atan2(-(Y - y_c), (X - x_c));
        theta = mod(theta, 2*pi);
        theta_ref = deg2rad(mod(RVrefAngleDeg(m), 360));

        ringType = lower(strtrim(slice_labels{m}));
        switch ringType
            case 'basal'
                nSeg = 6; sectorWidth = 2*pi/6; globalLabels = 1:6;
            case 'mid'
                nSeg = 6; sectorWidth = 2*pi/6; globalLabels = 7:12;
            case 'apical'
                nSeg = 4; sectorWidth = 2*pi/4; globalLabels = 13:16;
            otherwise
                error('Unknown slice label: %s', slice_labels{m});
        end

        if opts.clockwise
            angRel = mod(theta_ref - theta, 2*pi);
        else
            angRel = mod(theta - theta_ref, 2*pi);
        end

        rawSector = floor(mod(angRel + sectorWidth/2, 2*pi) / sectorWidth) + 1;
        rawSector(~Mask) = 0;

        localMap = zeros(size(Mask), 'like', rawSector);
        localMap(Mask) = rawSector(Mask);
        Mask_index_local(:,:,m) = localMap;

        % Final AHA remapping. This implements two things together:
        % 1) anticlockwise sector traversal in image-display coordinates
        % 2) a rotation of the starting reference sector so the resulting
        %    AHA identities align more closely with the supervisor's reference.
        switch ringType
            case 'basal'
                raw_to_global = [2 3 4 5 6 1];
            case 'mid'
                raw_to_global = [8 9 10 11 12 7];
            case 'apical'
                raw_to_global = [13 14 15 16];
        end

        globalMap = zeros(size(Mask), 'like', rawSector);
        for s = 1:nSeg
            globalMap(rawSector == s) = raw_to_global(s);
        end
        Mask_index_global(:,:,m) = globalMap;

        for s = 1:nSeg
            segID = raw_to_global(s);
            segMask = (globalMap == segID);
            validMask = segMask & isfinite(Img) & Img >= opts.validRange(1) & Img <= opts.validRange(2);
            raw_vals = double(Img(validMask));
            raw_vals = raw_vals(:);

            is_zero      = abs(raw_vals) <= opts.zero_tol;
            is_nan       = isnan(raw_vals);
            is_inf       = isinf(raw_vals);
            is_nonfinite = ~isfinite(raw_vals);
            is_bad       = is_zero | is_nan | is_inf | is_nonfinite;

            clean_vals = raw_vals(~is_bad);
            Segmentpix{segID,m} = clean_vals;

            n_total_raw   = numel(raw_vals);
            n_zero        = nnz(is_zero);
            n_nan         = nnz(is_nan);
            n_inf         = nnz(is_inf);
            n_nonfinite   = nnz(is_nonfinite);
            n_valid_clean = numel(clean_vals);
            n_bad         = nnz(is_bad);
            frac_bad = 0;
            if n_total_raw > 0, frac_bad = n_bad / n_total_raw; end

            qa_info.n_total_raw(segID,m)   = n_total_raw;
            qa_info.n_zero(segID,m)        = n_zero;
            qa_info.n_nan(segID,m)         = n_nan;
            qa_info.n_inf(segID,m)         = n_inf;
            qa_info.n_nonfinite(segID,m)   = n_nonfinite;
            qa_info.n_valid_clean(segID,m) = n_valid_clean;
            qa_info.frac_bad(segID,m)      = frac_bad;
            qa_info.flag_bad(segID,m)      = (n_bad >= opts.abs_bad_threshold) || (frac_bad >= opts.frac_bad_threshold);
        end
    end

    for segID = 1:16
        vals_all = [];
        for m = 1:nSlices
            vals_all = [vals_all; Segmentpix{segID,m}(:)]; %#ok<AGROW>
        end
        if isempty(vals_all)
            stats(:,segID) = [NaN; NaN; NaN; NaN; 0];
        else
            stats(:,segID) = [mean(abs(vals_all),'omitnan'); mean(vals_all,'omitnan'); std(vals_all,'omitnan'); median(vals_all,'omitnan'); numel(vals_all)];
        end
    end

    flagged_any = any(qa_info.flag_bad(:));
    if flagged_any
        report_file = fullfile(save_path, [filename sprintf('_segment_QA_phase_%02d.txt', ph)]);
        fid = fopen(report_file, 'w');
        if fid ~= -1
            fprintf(fid, 'AHA-corrected Segment QA report\n');
            fprintf(fid, 'File: %s\n', filename);
            fprintf(fid, 'Phase: %d\n', ph);
            fprintf(fid, 'Date: %s\n\n', datestr(now));
            for segID = 1:16
                for m = 1:nSlices
                    if qa_info.flag_bad(segID,m)
                        fprintf(fid, 'Segment %d | selected slice %d\n', segID, m);
                        fprintf(fid, '  total raw values   : %d\n', qa_info.n_total_raw(segID,m));
                        fprintf(fid, '  zero count         : %d\n', qa_info.n_zero(segID,m));
                        fprintf(fid, '  NaN count          : %d\n', qa_info.n_nan(segID,m));
                        fprintf(fid, '  Inf count          : %d\n', qa_info.n_inf(segID,m));
                        fprintf(fid, '  nonfinite count    : %d\n', qa_info.n_nonfinite(segID,m));
                        fprintf(fid, '  valid clean values : %d\n', qa_info.n_valid_clean(segID,m));
                        fprintf(fid, '  bad fraction       : %.4f\n\n', qa_info.frac_bad(segID,m));
                    end
                end
            end
            fclose(fid);
        else
            report_file = '';
        end
    end
end

function out = compute_bullseye_by_ring(metricMaps, segLabelMaps, baseIdx, midIdx, apexIdx, summaryFcn)
    if ~isequal(size(metricMaps), size(segLabelMaps))
        error('metricMaps and segLabelMaps must have identical size.');
    end
    targetSegs = 1:16;
    out.all16  = aggregate_over_slices(metricMaps, segLabelMaps, 1:size(metricMaps,3), targetSegs, summaryFcn);
    out.base16 = aggregate_over_slices(metricMaps, segLabelMaps, baseIdx, targetSegs, summaryFcn);
    out.mid16  = aggregate_over_slices(metricMaps, segLabelMaps, midIdx, targetSegs, summaryFcn);
    out.apex16 = aggregate_over_slices(metricMaps, segLabelMaps, apexIdx, targetSegs, summaryFcn);
    out.meta.baseSlices = baseIdx;
    out.meta.midSlices  = midIdx;
    out.meta.apexSlices = apexIdx;
    out.meta.summaryFcn = func2str(summaryFcn);
    out.meta.note = 'AHA-aligned 16-segment bull''s-eye from selected basal/mid/apical slices at one phase';
end

function segVals = aggregate_over_slices(metricMaps, segLabelMaps, sliceIdx, targetSegs, summaryFcn)
    segVals = nan(numel(targetSegs),1);
    sliceIdx = sliceIdx(sliceIdx >= 1 & sliceIdx <= size(metricMaps,3));
    if isempty(sliceIdx), return; end
    for i = 1:numel(targetSegs)
        seg = targetSegs(i);
        vals = [];
        for s = sliceIdx
            m = metricMaps(:,:,s);
            lab = segLabelMaps(:,:,s);
            mask = (lab == seg) & isfinite(m);
            if any(mask(:))
                vals = [vals; m(mask)]; %#ok<AGROW>
            end
        end
        if ~isempty(vals)
            segVals(i) = summaryFcn(vals);
        end
    end
end

function plot_bullseye_aha16_AHA(v16, titleStr)
    v = v16(:);
    if numel(v) ~= 16, error('Input must be a 16-element vector.'); end

    cla; hold on; axis equal off;
    title(titleStr, 'Interpreter', 'none');
    finiteVals = v(isfinite(v));
    if isempty(finiteVals)
        clim = [0 1];
    else
        clim = [min(finiteVals) max(finiteVals)];
        if clim(1) == clim(2)
            pad = max(abs(clim(1))*0.05, 1e-6);
            clim = clim + [-pad pad];
        end
    end
    caxis(clim); colormap(parula); colorbar;

    rOuter = 3.0; rMid2 = 2.0; rMid1 = 1.0; rInner = 0.35;

    draw_ring_segments(v(1:6),  1:6,   rMid2, rOuter, 6);
    draw_ring_segments(v(7:12), 7:12,  rMid1, rMid2, 6);
    draw_ring_segments(v(13:16),13:16, rInner, rMid1, 4);

    th = linspace(0, 2*pi, 200);
    patch(rInner*cos(th), rInner*sin(th), 'w', 'EdgeColor', 'k', 'LineWidth', 1);
    xlim([-3.5 4.2]); ylim([-3.5 3.5]); hold off;
end

function draw_ring_segments(vals, segNums, rIn, rOut, nSeg)
    sector = 2*pi / nSeg;
    for k = 1:nSeg
        segValue = vals(k);
        segNum   = segNums(k);
        thetaCenter = pi/2 - (k-1)*sector;  % top then clockwise
        theta1 = thetaCenter + sector/2;
        theta2 = thetaCenter - sector/2;
        h = draw_sector(rIn, rOut, theta1, theta2, segValue);
        if isnan(segValue)
            set(h, 'FaceColor', [0.85 0.85 0.85]);
        end
        rText = 0.5*(rIn + rOut);
        xText = rText*cos(thetaCenter);
        yText = rText*sin(thetaCenter);
        if isnan(segValue)
            txt = sprintf('%d\nNaN', segNum);
        else
            txt = sprintf('%d\n%.3g', segNum, segValue);
        end
        text(xText, yText, txt, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',9, 'FontWeight','bold');
    end
end

function h = draw_sector(rIn, rOut, theta1, theta2, cval)
    th = linspace(theta1, theta2, 80);
    x = [rOut*cos(th), fliplr(rIn*cos(th))];
    y = [rOut*sin(th), fliplr(rIn*sin(th))];
    if isnan(cval)
        h = patch(x, y, [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1);
    else
        h = patch(x, y, cval, 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'flat');
    end
end

function overlay_aha_labels(ax, maskIndexSlice, sliceType)
    labels = unique(maskIndexSlice(:));
    labels(labels==0) = [];
    switch lower(sliceType)
        case 'basal'
            mapNames = {'1 Ant','2 AS','3 IS','4 Inf','5 IL','6 AL'};
            labelKeys = 1:6;
        case 'mid'
            mapNames = {'7 Ant','8 AS','9 IS','10 Inf','11 IL','12 AL'};
            labelKeys = 7:12;
        case 'apical'
            mapNames = {'13 Ant','14 Sept','15 Inf','16 Lat'};
            labelKeys = 13:16;
        otherwise
            mapNames = {}; labelKeys = [];
    end
    nameMap = containers.Map(num2cell(labelKeys), mapNames);
    for k = labels(:)'
        BW = maskIndexSlice==k;
        if any(BW(:))
            S = regionprops(BW,'Centroid','Area');
            [~,iMax] = max([S.Area]);
            c = S(iMax).Centroid;
            if isKey(nameMap,k), txt = nameMap(k); else, txt = sprintf('%d',k); end
            text(ax, c(1), c(2), txt, 'Color','k', 'FontSize',10, 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end
end

function C = largest_component_props(BW)
    C = struct('Centroid',[NaN NaN]);
    if ~any(BW(:)), return; end
    CC = bwconncomp(BW);
    if CC.NumObjects==0, return; end
    stats = regionprops(CC,'Area','Centroid');
    [~,iMax] = max([stats.Area]);
    C = stats(iMax);
end

function [BW, ok, endo, epi] = get_or_edit_myo_roi_nozoom_click2edit_safe(imgMag, imgT2s, seedEndo, seedEpi, maskRef)
if nargin < 3, seedEndo = []; end
if nargin < 4, seedEpi = []; end
if nargin < 5, maskRef = []; end

fig = figure('Name','AHA ROI: Endocardial & Epicardial Sync','NumberTitle','off', ...
    'Units','normalized','Position',[0.05 0.1 0.9 0.75]);
tl  = tiledlayout(fig,1,2,'TileSpacing','compact');

ax1 = nexttile(tl,1);
imshow(imgMag,[],'Parent',ax1); hold(ax1,'on');
title(ax1, '1. Draw/Adjust BLUE (Endo) then double-click');
if ~isempty(maskRef) && any(maskRef(:))
    visboundaries(ax1, maskRef, 'Color', [1 1 0], 'LineWidth', 0.5);
end

ax2 = nexttile(tl,2);
imagesc(ax2, imgT2s); axis(ax2,'image'); axis(ax2,'off'); colormap(ax2,'jet'); caxis(ax2,[0 80]); colorbar(ax2); hold(ax2,'on');
title(ax2,'T2* Live Mirror');
if ~isempty(maskRef) && any(maskRef(:))
    visboundaries(ax2, maskRef, 'Color', [1 1 0], 'LineWidth', 0.5);
end
linkaxes([ax1 ax2],'xy');

[hEndo, ~] = create_synced_poly(ax1, ax2, seedEndo, [0 1 1]);
wait(hEndo);
title(ax1, '2. Draw/Adjust PURPLE (Epi) then double-click');
[hEpi, ~] = create_synced_poly(ax1, ax2, seedEpi, [1 0 1]);
wait(hEpi);

endo = hEndo.Position;
epi  = hEpi.Position;
[H,W] = size(imgMag);
BWendo = poly2mask(endo(:,1), endo(:,2), H, W);
BWepi  = poly2mask(epi(:,1),  epi(:,2),  H, W);
BW = logical(BWepi - BWendo > 0);
ok = ~isempty(endo) && ~isempty(epi);
if ishghandle(fig), close(fig); end
end

function [hROI, hMirror] = create_synced_poly(axSrc, axTarget, seed, color)
if ~isempty(seed)
    hROI = drawpolygon(axSrc,'Position',seed,'Color',color,'FaceAlpha',0.1);
else
    hROI = drawpolygon(axSrc,'Color',color,'FaceAlpha',0.1);
end
hMirror = plot(axTarget, NaN, NaN, 'Color', color, 'LineWidth', 1.5);
addlistener(hROI, 'MovingROI', @(src, evt) update_mirror_line(src, hMirror));
update_mirror_line(hROI, hMirror);
end

function update_mirror_line(roiObj, lineObj)
p = roiObj.Position;
if ~isempty(p)
    closedP = [p; p(1,:)];
    set(lineObj, 'XData', closedP(:,1), 'YData', closedP(:,2));
end
end

function [pt, accepted] = get_or_edit_point(imgGray, seedPt)
if nargin<2 || isempty(seedPt), seedPt = [NaN NaN]; end
fig = figure('Name','Anterior RV insertion point','NumberTitle','off','Units','normalized','Position',[0.33 0.33 0.34 0.34]);
ax = axes(fig); imshow(imgGray,[],'Parent',ax); hold(ax,'on');
zoom(fig,'off'); pan(fig,'off');
if ~any(isnan(seedPt))
    plot(ax, seedPt(1),seedPt(2),'r+','MarkerSize',10,'LineWidth',1);
    title(ax,'Click the ANTERIOR RV insertion point, or press A to accept seed; Q to cancel');
else
    title(ax,'Click the ANTERIOR RV insertion point; Q to cancel');
end
accepted = false; pt = seedPt;
set(fig,'WindowButtonDownFcn',@onClick); set(fig,'KeyPressFcn',@onKey); uiwait(fig);
    function onClick(~,~)
        cp = get(ax,'CurrentPoint'); pt = cp(1,1:2); accepted = true; uiresume(fig); close(fig);
    end
    function onKey(~,evt)
        switch lower(evt.Key)
            case 'a'
                if ~any(isnan(seedPt)), accepted = true; uiresume(fig); close(fig); end
            case {'q','escape'}
                accepted = false; uiresume(fig); close(fig);
        end
    end
end

function save_current_figure(fname)
if exist('exportgraphics','file') == 2
    exportgraphics(gcf, fname, 'Resolution', 300);
else
    saveas(gcf, fname);
end
end
