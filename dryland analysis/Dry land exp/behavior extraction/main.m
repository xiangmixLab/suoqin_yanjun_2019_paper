%dry land behav mean
%orilocation: data location
%conditionName: name for each orilocation
%behavVidprefix: the first common part of all behav videos
%mrange: which orilocation need be processed
%ROImethod: 
%'autoROI' - premitive automatic ROI detection, not always work;
%'manualROI' - manual ROI selection
%'ROIlist' - manually define all ROI prior processing
%ROIlist: when using 'ROIlist' option, input the predefined ROIs
%objInput: 
%'Y': select Objects (for OLM or ORM)
%'List': pre-select all objects prior processing
%objList: when using 'List' option, input the predefined objs
[behavname]=Miniscope_behav_extraction_auto_dry_land(orilocation,conditionName,behavVidprefix,mrange,ROImethod,ROIlist,objInput,objlist)
