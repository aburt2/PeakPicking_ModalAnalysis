classdef modalAnalyser < handle
    %Peack Picking algorithm for modal expansion
    %   Detailed explanation goes here

    properties (Access = public)
        modalFrequency
        modalAmplitudes
        poleIdx
        impMeasurement
        v
        x
        M
        d
        damping
        fittingBandwidth
        simImpedance
    end

    methods
        function obj = modalAnalyser(damping,fittingBandwidth,impMeasurement)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.modalFrequency = [];
            obj.v = [];
            obj.modalAmplitudes = [];
            obj.M = [];
            obj.d = [];
            obj.simImpedance = [];
            obj.damping = damping;
            obj.fittingBandwidth = fittingBandwidth;
            obj.impMeasurement = impMeasurement;
        end

        function obj = addMode(obj,freqGuess)
            %ADDMODE Add new mode to modal expansion of the input frequency

            % Update frequency range
            obj.updateFrequencyRange(freqGuess)

            % Estimate modal frequency of new pole using frequency guess
            modeGuess = [freqGuess obj.damping];
            obj.modalFrequency = [obj.modalFrequency; modeGuess];
            obj.poleIdx = size(obj.modalFrequency,1); %location to add pole to in basis vector

            % Estimate frequency and damping of mode
            newMode = fminsearch(@(mode) obj.poleCostFunction(mode),modeGuess);
            obj.modalFrequency(obj.poleIdx,:) = newMode;
        end

        function obj = removeMode(obj,poleIndex)
            %REMOVEMODE Remove a mode from the basis

            % Remove mode from array
            obj.modalFrequency(poleIndex,:) = [];
            obj.modalAmplitudes(poleIndex,:) = [];

            % Update M matrix
            obj.updateMMatrix();

            % Make d matrixs
            obj.updateDMatrix();

            % Update amplitudes
            obj.modalAmplitudes = lsqminnorm(obj.M,obj.d);
        end
   
        function obj = refineIdentification(obj)
            %REFINEIDENTIFICATION Resolve equation 17 for each pole in the
            %basis

            % Sort poles in ascending order
            obj.modalFrequency = sort(obj.modalFrequency,1);

            % Start loop
            for n = 1:size(obj.modalFrequency,1)
                modeGuess = obj.modalFrequency(n,:);
                obj.poleIdx = n;

                % Estimate frequency and damping of mode
                newMode = fminsearch(@(mode) obj.poleCostFunction(mode),modeGuess);
                obj.modalFrequency(obj.poleIdx,:) = newMode;
            end
        end
   
        function obj = changeBandwidth(obj,newBandwidth)
            %CHANGEBANDWIDTH Update the fitting bandwidth with new value
            %   Detailed explanation goes here
            obj.fittingBandwidth = newBandwidth;

            % Recompute v
            for n = 1:size(obj.modalFrequency,1)
                wn = obj.modalFrequency(n,1);
                obj.updateFrequencyRange(wn)
            end
        end
    
        function obj = updateMMatrix(obj)
            %UPDATEMMATRIX update M matrix based on new mode guess

            % Get the frequencies and damping
            w = obj.modalFrequency(:,1);
            zeta = obj.modalFrequency(:,2);
            
            % Make empty matrix
            Ma = zeros(length(obj.v),length(w));
            Mb = zeros(length(obj.v),length(w));
            
            % Common denominator
            for n = 1:length(w)
                % Get frequency and damping
                wn = w(n);
                zetan = zeta(n);

                % Compute common denominator
                denom = wn.^2-obj.v.^2+1j.*2.*zetan.*wn.*obj.v;
                
                % Form Ma and Mb matrices
                Ma(:,n) = 1j.*obj.v./denom;
                Mb(:,n) = 1./denom;
            end
                
            % Create M matrix
            obj.M = [real(Ma) real(Mb);imag(Ma) imag(Mb)];
        end

        function obj = updateDMatrix(obj)
            %IMPEDANCEPLOT get an array of frequencies and associated
            %magnitudes given the current poles
            
            % Get Frequency measurements
            logicIdx = ismember(obj.impMeasurement,obj.v);

            % Get magnitudes
            temp = obj.impMeasurement(logicIdx,2);
            
            % Save array
            obj.d = [real(temp); imag(temp)];
        end

        function obj = impedancePlot(obj, freq)
            %IMPEDANCEPLOT get an array of frequencies and associated
            %magnitudes given the current poles
            
            % Get modes and damping
            w = obj.modalFrequency(:,1);
            zeta = obj.modalFrequency(:,2);

            % Split complex and real amplitudes
            k = length(obj.modalAmplitudes);
            A = obj.modalAmplitudes(1:k/2);
            B = obj.modalAmplitudes(k/2+1:end);

            % Create empty vectors
            temp = zeros(length(freq),2);
            temp(:,1) = real(freq);

            % Compute magnitudes
            for n = 1:length(w)
                %compute numerator
                num = 1j*freq*A(n)+B(n);
                denom = w(n)^2 - freq.^2 + 1j*2*zeta(n).*freq*w(n);

                temp(:,2) = temp(:,2) + num./denom;
            end

            % Save array
            obj.simImpedance = temp;
        end

        function obj = updateFrequencyRange(obj,freqGuess)
            %UPDATEFREQUENCYRANGE Update the frequency range v with new
            %guess
            % Get just the frequencies
            freqMeasurement = obj.impMeasurement(:,1);

            % Find the minimum and maximum ranges
            freqRangeMin = freqGuess - obj.fittingBandwidth/2;
            freqRangeMax = freqGuess + obj.fittingBandwidth/2;

            freqRange = freqMeasurement(and(freqMeasurement>freqRangeMin,freqMeasurement<freqRangeMax));

            obj.v = union(obj.v,freqRange);
        end

        function cost = poleCostFunction(obj,modeGuess)
            %POLECOSTFUNCTION Cost function for estimating new pole

            % Update modal frequencies with guess
            obj.modalFrequency(obj.poleIdx,:) = modeGuess;

            % Get matrices
            obj.updateMMatrix();

            % Make d matrixs
            obj.updateDMatrix();
            
            % Guess amplitudes
            obj.modalAmplitudes = lsqminnorm(obj.M,obj.d);

            % return value
            cost = norm(obj.M*obj.modalAmplitudes - obj.d);
        end
    end
end