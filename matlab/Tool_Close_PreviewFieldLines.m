function out=Tool_Close_PreviewFieldLines(h)
for i=1:length(h)
    for j = 1:length(h{i})
        set(h{i}(j),'visible','off');
    end
end