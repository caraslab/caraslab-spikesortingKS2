function dickey = computeDailyDickey(spiketimes)
warning('off','stats:gmdistribution:FailedToConverge');
dickey = cell(size(spiketimes));

fprintf('Computing daily dickey\n');
for iid=1:length(spiketimes)
    fprintf('%d',iid);
    dickey{iid} = zeros(8,length(spiketimes{iid}));
    
    for iin=1:length(spiketimes{iid})
        fprintf('.');
        isi = log(diff(spiketimes{iid}{iin}));
        try 
            % This is an odd way to fit this distribution, given that we
            % don't actually believe the ISIS are log normal.  It would be
            % more logical to use nonlinear least squares curve fitting but
            % I have to do it this way to be the same as Dickey
            %
            % This is not very many max iterations but I have checked it
            % and they fit just fine.  Also the scores that are used for
            % classification look if anything better than the plot in
            % dickey's paper
            model = gmdistribution.fit(isi,3,'Start',struct('mu',[-6 -3.5 0]','Sigma',shiftdim([.5 .9 10],-1),'PComponents',[.02 .6 .38]'),'options',struct('MaxIter',10));
            dickey{iid}(:,iin) = [model.mu(:)' model.Sigma(:)' model.PComponents(1:2)];
        catch
        end
        
%         edges = min(isi):range(isi)/100:max(isi);
%         edgeCenters = (edges(1:end-1)+edges(2:end))/2;
%         clf, hold on
%         h = histc(isi,edges);
%         h = h(1:end-1);
%         h = h./sum(h)./diff(edges)';
%         bar(edgeCenters,h);
%         plot(edgeCenters,pdf(model,edgeCenters'),'r','LineWidth',2);
%         model = gmdistribution.fit(isi,3,'Start',struct('mu',[-6 -3.5 0]','Sigma',shiftdim([.5 .9 10],-1),'PComponents',[.02 .6 .38]'),'options',struct('MaxIter',100));
%         plot(edgeCenters,pdf(model,edgeCenters'),'g','LineWidth',2);
%         keyboard;
    end
    fprintf('\n');
end

warning('on','stats:gmdistribution:FailedToConverge');
end