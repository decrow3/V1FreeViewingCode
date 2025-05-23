function cmap = getcolormap(subject, plotit)
% cmap = getcolormap(subject)

if nargin < 2
    plotit = true;
end

if nargin < 1
    subject = 'allen';
end

switch subject
    case 'gru'
        cmap = [[217, 223, 228]; ...
            [182, 202, 215]; ...
            [144, 181, 207]; ...
            [105, 162, 201]; ...
            [66, 143, 195]; ...
            [51, 121, 169]; ...
            [44, 99, 136]; ...
            [39, 77, 103]; ...
            [34, 57, 72]; ...
            [26, 37, 44]]/255;

    case 'brie'

        cmap = [[209, 225, 213]; ...
            [166, 213, 178]; ...
            [119, 208, 141]; ...
            [70, 207, 104]; ...
            [40, 190, 78]; ...
            [31, 156, 63]; ...
            [27, 122, 51]; ...
            [25, 90, 41]; ...
            [23, 61, 32]; ...
            [19, 36, 23]]/255;

    case 'allen'

        cmap = [[233, 226, 221]; ...
            [224, 203, 189]; ...
            [220, 180, 152]; ...
            [218, 156, 112]; ...
            [216, 133, 73]; ...
            [206, 110, 41]; ...
            [169, 93, 38]; ...
            [129, 75, 36]; ...
            [91, 58, 34]; ...
            [55, 40, 29]]/255;

end

if plotit
    figure; clf
    n = size(cmap,1);
    x = linspace(0, 10, 100);
    for i = 1:n
        plot(x, i*x, 'Color', cmap(i,:), 'Linewidth', 5); hold on
    end
    title(subject)
end
