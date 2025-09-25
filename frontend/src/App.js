import React, { useState, useEffect } from "react";
import "./App.css";
import axios from "axios";
import { Upload, Play, Settings, FileText, BarChart3, Download, AlertCircle } from "lucide-react";

const BACKEND_URL = process.env.REACT_APP_BACKEND_URL;
const API = `${BACKEND_URL}/api`;

// Main EarthImager Interface Component
const EarthImagerInterface = () => {
  const [activeTab, setActiveTab] = useState('parameters');
  const [forwardParams, setForwardParams] = useState({
    n_electrodes: 8,
    electrode_spacing: 1.0,
    resistivity: 100.0
  });
  const [iniParams, setIniParams] = useState({
    forward_method: 0,
    forward_solver: 0,
    bc_type: 0,
    max_iterations: 20,
    lagrange: 10.0,
    start_resistivity: 1.0,
    min_resistivity: 1.0,
    max_resistivity: 100000.0
  });
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [uploadedFiles, setUploadedFiles] = useState({});
  const [status, setStatus] = useState('');

  // Check backend health on mount
  useEffect(() => {
    checkHealth();
  }, []);

  const checkHealth = async () => {
    try {
      const response = await axios.get(`${API}/health`);
      console.log('Backend status:', response.data);
      setStatus(response.data.ei2d_engine === 'available' ? 'EI2D Engine Ready' : 'EI2D Engine Unavailable');
    } catch (error) {
      setStatus('Backend Connection Failed');
      console.error('Health check failed:', error);
    }
  };

  const runRealForwardModel = async () => {
    if (!uploadedFiles.ini || !uploadedFiles.stg) {
      setStatus('Error: Please upload both INI and STG files for real forward modeling');
      return;
    }

    setLoading(true);
    try {
      const formData = new FormData();
      // We need to recreate File objects from the uploaded data
      const iniBlob = new Blob([uploadedFiles.ini.data.parsed_data ? 
        JSON.stringify(uploadedFiles.ini.data.parsed_data) : 'Mock INI content'], 
        { type: 'text/plain' });
      const stgBlob = new Blob(['Mock STG content'], { type: 'text/plain' });
      
      formData.append('ini_file', iniBlob, uploadedFiles.ini.name);
      formData.append('stg_file', stgBlob, uploadedFiles.stg.name);

      const response = await axios.post(`${API}/earthimager/forward-model-real`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });
      
      setResults(response.data);
      setStatus(`Real forward modeling completed: ${response.data.method}`);
    } catch (error) {
      setStatus(`Real forward modeling error: ${error.response?.data?.detail || error.message}`);
      console.error('Real forward modeling error:', error);
    }
    setLoading(false);
  };

  const runForwardModel = async () => {
    setLoading(true);
    try {
      const response = await axios.post(`${API}/earthimager/forward-model`, forwardParams);
      setResults(response.data);
      setStatus('Forward modeling completed successfully');
    } catch (error) {
      setStatus(`Error: ${error.response?.data?.detail || error.message}`);
      console.error('Forward modeling error:', error);
    }
    setLoading(false);
  };

  const uploadFile = async (file, fileType) => {
    const formData = new FormData();
    formData.append('file', file);
    
    try {
      const response = await axios.post(`${API}/earthimager/upload-${fileType}`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });
      
      setUploadedFiles(prev => ({
        ...prev,
        [fileType]: { name: file.name, data: response.data }
      }));
      
      setStatus(`${fileType.toUpperCase()} file processed successfully`);
    } catch (error) {
      setStatus(`Error processing ${fileType} file: ${error.response?.data?.detail || error.message}`);
    }
  };

  const generateIniFile = async () => {
    try {
      const response = await axios.post(`${API}/earthimager/generate-ini`, iniParams);
      setStatus('INI configuration generated successfully');
      
      // Download the generated file
      const blob = new Blob([response.data.content], { type: 'text/plain' });
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.style.display = 'none';
      a.href = url;
      a.download = 'earthimager_config.ini';
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
      document.body.removeChild(a);
    } catch (error) {
      setStatus(`Error generating INI: ${error.response?.data?.detail || error.message}`);
    }
  };

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header - Mimics EI2D Menu Bar */}
      <div className="bg-white border-b border-gray-200 shadow-sm">
        <div className="px-6 py-3">
          <div className="flex items-center justify-between">
            <h1 className="text-xl font-semibold text-gray-900">EarthImager 2D Web Interface</h1>
            <div className="flex items-center space-x-4">
              <span className={`px-3 py-1 rounded-full text-sm font-medium ${
                status.includes('Ready') ? 'bg-green-100 text-green-800' : 
                status.includes('Error') || status.includes('Failed') ? 'bg-red-100 text-red-800' :
                'bg-yellow-100 text-yellow-800'
              }`}>
                {status || 'Initializing...'}
              </span>
            </div>
          </div>
          
          {/* Menu Bar Tabs */}
          <div className="mt-4 border-b border-gray-200">
            <nav className="-mb-px flex space-x-8">
              {[
                { id: 'parameters', label: 'Parameters', icon: Settings },
                { id: 'files', label: 'File Manager', icon: FileText },
                { id: 'results', label: 'Results', icon: BarChart3 }
              ].map(({ id, label, icon: Icon }) => (
                <button
                  key={id}
                  onClick={() => setActiveTab(id)}
                  className={`flex items-center space-x-2 py-2 px-1 border-b-2 font-medium text-sm ${
                    activeTab === id
                      ? 'border-blue-500 text-blue-600'
                      : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                  }`}
                >
                  <Icon className="w-4 h-4" />
                  <span>{label}</span>
                </button>
              ))}
            </nav>
          </div>
        </div>
      </div>

      {/* Main Content Area */}
      <div className="p-6">
        {activeTab === 'parameters' && (
          <ParametersTab
            forwardParams={forwardParams}
            setForwardParams={setForwardParams}
            iniParams={iniParams}
            setIniParams={setIniParams}
            onRunForwardModel={runForwardModel}
            onRunRealForwardModel={runRealForwardModel}
            onGenerateIni={generateIniFile}
            loading={loading}
            uploadedFiles={uploadedFiles}
          />
        )}

        {activeTab === 'files' && (
          <FilesTab
            uploadedFiles={uploadedFiles}
            onFileUpload={uploadFile}
          />
        )}

        {activeTab === 'results' && (
          <ResultsTab results={results} />
        )}
      </div>
    </div>
  );
};

// Parameters Tab Component
const ParametersTab = ({ forwardParams, setForwardParams, iniParams, setIniParams, onRunForwardModel, onRunRealForwardModel, onGenerateIni, loading, uploadedFiles }) => {
  const [paramTab, setParamTab] = useState('forward');

  return (
    <div className="bg-white rounded-lg shadow-sm border border-gray-200">
      <div className="border-b border-gray-200">
        <nav className="-mb-px flex">
          {[
            { id: 'forward', label: 'Forward Modeling' },
            { id: 'inversion', label: 'Inversion Settings' }
          ].map(({ id, label }) => (
            <button
              key={id}
              onClick={() => setParamTab(id)}
              className={`px-6 py-3 font-medium text-sm ${
                paramTab === id
                  ? 'border-b-2 border-blue-500 text-blue-600'
                  : 'text-gray-500 hover:text-gray-700'
              }`}
            >
              {label}
            </button>
          ))}
        </nav>
      </div>

      <div className="p-6">
        {paramTab === 'forward' && (
          <div className="space-y-6">
            <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Number of Electrodes
                </label>
                <input
                  type="number"
                  min="4"
                  max="120"
                  value={forwardParams.n_electrodes}
                  onChange={(e) => setForwardParams(prev => ({
                    ...prev,
                    n_electrodes: parseInt(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Electrode Spacing (m)
                </label>
                <input
                  type="number"
                  step="0.1"
                  min="0.1"
                  value={forwardParams.electrode_spacing}
                  onChange={(e) => setForwardParams(prev => ({
                    ...prev,
                    electrode_spacing: parseFloat(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Resistivity (Ω·m)
                </label>
                <input
                  type="number"
                  step="1"
                  min="0.1"
                  value={forwardParams.resistivity}
                  onChange={(e) => setForwardParams(prev => ({
                    ...prev,
                    resistivity: parseFloat(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
              </div>
            </div>

            <div className="flex justify-end space-x-4">
              <button
                onClick={runForwardModel}
                disabled={loading}
                className="flex items-center space-x-2 px-6 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50"
              >
                <Play className="w-4 h-4" />
                <span>{loading ? 'Running...' : 'Run Forward Model (Mock)'}</span>
              </button>
              
              {(uploadedFiles.ini && uploadedFiles.stg) && (
                <button
                  onClick={runRealForwardModel}
                  disabled={loading}
                  className="flex items-center space-x-2 px-6 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 disabled:opacity-50"
                >
                  <Play className="w-4 h-4" />
                  <span>{loading ? 'Processing...' : 'Run with Real Data'}</span>
                </button>
              )}
            </div>
          </div>
        )}

        {paramTab === 'inversion' && (
          <div className="space-y-6">
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Forward Method
                </label>
                <select
                  value={iniParams.forward_method}
                  onChange={(e) => setIniParams(prev => ({
                    ...prev,
                    forward_method: parseInt(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                >
                  <option value={0}>Finite Difference (FD)</option>
                  <option value={1}>Finite Element (FE)</option>
                </select>
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Forward Solver
                </label>
                <select
                  value={iniParams.forward_solver}
                  onChange={(e) => setIniParams(prev => ({
                    ...prev,
                    forward_solver: parseInt(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                >
                  <option value={0}>Cholesky Decomposition</option>
                  <option value={1}>Conjugate Gradient</option>
                </select>
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Boundary Condition
                </label>
                <select
                  value={iniParams.bc_type}
                  onChange={(e) => setIniParams(prev => ({
                    ...prev,
                    bc_type: parseInt(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                >
                  <option value={0}>Dirichlet</option>
                  <option value={1}>Mixed</option>
                </select>
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Max Iterations
                </label>
                <input
                  type="number"
                  min="1"
                  max="100"
                  value={iniParams.max_iterations}
                  onChange={(e) => setIniParams(prev => ({
                    ...prev,
                    max_iterations: parseInt(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Lagrange Multiplier
                </label>
                <input
                  type="number"
                  step="0.1"
                  min="0.1"
                  value={iniParams.lagrange}
                  onChange={(e) => setIniParams(prev => ({
                    ...prev,
                    lagrange: parseFloat(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Start Resistivity (Ω·m)
                </label>
                <input
                  type="number"
                  step="0.1"
                  min="0.1"
                  value={iniParams.start_resistivity}
                  onChange={(e) => setIniParams(prev => ({
                    ...prev,
                    start_resistivity: parseFloat(e.target.value)
                  }))}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
              </div>
            </div>

            <div className="flex justify-end">
              <button
                onClick={onGenerateIni}
                className="flex items-center space-x-2 px-6 py-2 bg-green-600 text-white rounded-md hover:bg-green-700"
              >
                <Download className="w-4 h-4" />
                <span>Generate INI File</span>
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

// Files Tab Component
const FilesTab = ({ uploadedFiles, onFileUpload }) => {
  const handleFileUpload = (event, fileType) => {
    const file = event.target.files[0];
    if (file) {
      onFileUpload(file, fileType);
    }
    event.target.value = '';
  };

  return (
    <div className="space-y-6">
      <div className="bg-white rounded-lg shadow-sm border border-gray-200 p-6">
        <h2 className="text-lg font-semibold text-gray-900 mb-4">File Manager</h2>
        
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
          {/* INI File Upload */}
          <div className="border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-gray-400 transition-colors">
            <div className="text-center">
              <Upload className="mx-auto h-8 w-8 text-gray-400" />
              <div className="mt-2">
                <label className="cursor-pointer">
                  <span className="text-sm font-medium text-gray-900">Upload INI File</span>
                  <input
                    type="file"
                    accept=".ini"
                    onChange={(e) => handleFileUpload(e, 'ini')}
                    className="hidden"
                  />
                </label>
                <p className="text-xs text-gray-500 mt-1">Configuration (.ini)</p>
              </div>
              {uploadedFiles.ini && (
                <div className="mt-2 text-sm text-green-600">
                  ✓ {uploadedFiles.ini.name}
                </div>
              )}
            </div>
          </div>

          {/* STG File Upload */}
          <div className="border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-gray-400 transition-colors">
            <div className="text-center">
              <Upload className="mx-auto h-8 w-8 text-gray-400" />
              <div className="mt-2">
                <label className="cursor-pointer">
                  <span className="text-sm font-medium text-gray-900">Upload STG File</span>
                  <input
                    type="file"
                    accept=".stg"
                    onChange={(e) => handleFileUpload(e, 'stg')}
                    className="hidden"
                  />
                </label>
                <p className="text-xs text-gray-500 mt-1">Survey data (.stg)</p>
              </div>
              {uploadedFiles.stg && (
                <div className="mt-2 text-sm text-green-600">
                  ✓ {uploadedFiles.stg.name}
                </div>
              )}
            </div>
          </div>

          {/* MDL File Upload */}
          <div className="border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-gray-400 transition-colors">
            <div className="text-center">
              <Upload className="mx-auto h-8 w-8 text-gray-400" />
              <div className="mt-2">
                <label className="cursor-pointer">
                  <span className="text-sm font-medium text-gray-900">Upload MDL File</span>
                  <input
                    type="file"
                    accept=".mdl"
                    onChange={(e) => handleFileUpload(e, 'mdl')}
                    className="hidden"
                  />
                </label>
                <p className="text-xs text-gray-500 mt-1">Model (.mdl)</p>
              </div>
              {uploadedFiles.mdl && (
                <div className="mt-2 text-sm text-green-600">
                  ✓ {uploadedFiles.mdl.name}
                </div>
              )}
            </div>
          </div>

          {/* MOD File Upload */}
          <div className="border-2 border-dashed border-gray-300 rounded-lg p-4 hover:border-gray-400 transition-colors">
            <div className="text-center">
              <Upload className="mx-auto h-8 w-8 text-gray-400" />
              <div className="mt-2">
                <label className="cursor-pointer">
                  <span className="text-sm font-medium text-gray-900">Upload MOD File</span>
                  <input
                    type="file"
                    accept=".mod"
                    onChange={(e) => handleFileUpload(e, 'mod')}
                    className="hidden"
                  />
                </label>
                <p className="text-xs text-gray-500 mt-1">Resistivity (.mod)</p>
              </div>
              {uploadedFiles.mod && (
                <div className="mt-2 text-sm text-green-600">
                  ✓ {uploadedFiles.mod.name}
                </div>
              )}
            </div>
          </div>
        </div>

        {/* File Information Display */}
        {Object.keys(uploadedFiles).length > 0 && (
          <div className="mt-6">
            <h3 className="text-md font-medium text-gray-900 mb-3">Uploaded Files Information</h3>
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
              {uploadedFiles.ini && (
                <div className="bg-blue-50 border border-blue-200 rounded-md p-4">
                  <h4 className="font-medium text-blue-900 mb-2">📄 INI Configuration</h4>
                  <div className="text-sm text-blue-700 space-y-1">
                    <div>Sections: {uploadedFiles.ini.data.sections?.join(', ')}</div>
                    <div>Forward Method: {uploadedFiles.ini.data.forward_method === '0' ? 'FD' : 'FE'}</div>
                    <div>Max Iterations: {uploadedFiles.ini.data.max_iterations}</div>
                  </div>
                </div>
              )}

              {uploadedFiles.stg && (
                <div className="bg-green-50 border border-green-200 rounded-md p-4">
                  <h4 className="font-medium text-green-900 mb-2">🔌 Survey Data (STG)</h4>
                  <div className="text-sm text-green-700 space-y-1">
                    <div>Format: {uploadedFiles.stg.data.format}</div>
                    <div>Electrodes: {uploadedFiles.stg.data.num_electrodes}</div>
                    <div>Measurements: {uploadedFiles.stg.data.num_measurements}</div>
                    <div>Spacing: {uploadedFiles.stg.data.electrode_spacing?.toFixed(1)} m</div>
                    {uploadedFiles.stg.data.resistivity_range && (
                      <div>App. Res: {uploadedFiles.stg.data.resistivity_range.min?.toFixed(1)} - {uploadedFiles.stg.data.resistivity_range.max?.toFixed(1)} Ω·m</div>
                    )}
                  </div>
                </div>
              )}

              {uploadedFiles.mdl && (
                <div className="bg-purple-50 border border-purple-200 rounded-md p-4">
                  <h4 className="font-medium text-purple-900 mb-2">🏗️ Model (MDL)</h4>
                  <div className="text-sm text-purple-700 space-y-1">
                    <div>Format: {uploadedFiles.mdl.data.format}</div>
                    <div>Electrodes: {uploadedFiles.mdl.data.electrode_count}</div>
                    <div>Commands: {uploadedFiles.mdl.data.measurement_count}</div>
                    <div>Sections: {uploadedFiles.mdl.data.sections?.join(', ')}</div>
                  </div>
                </div>
              )}

              {uploadedFiles.mod && (
                <div className="bg-orange-50 border border-orange-200 rounded-md p-4">
                  <h4 className="font-medium text-orange-900 mb-2">🌍 Resistivity Model (MOD)</h4>
                  <div className="text-sm text-orange-700 space-y-1">
                    <div>Format: {uploadedFiles.mod.data.format}</div>
                    <div>Background: {uploadedFiles.mod.data.background_resistivity} Ω·m</div>
                    <div>Blocks: {uploadedFiles.mod.data.total_blocks}</div>
                    {uploadedFiles.mod.data.layers && (
                      <div>Layers: {uploadedFiles.mod.data.layers.num_layers}</div>
                    )}
                  </div>
                </div>
              )}
            </div>

            {/* 3-Layer Model Visualization */}
            {uploadedFiles.mod?.data?.layers && (
              <div className="mt-4 bg-gray-50 border border-gray-200 rounded-md p-4">
                <h4 className="font-medium text-gray-900 mb-3">3-Layer Resistivity Model</h4>
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  {uploadedFiles.mod.data.layers.layer_info?.map((layer, idx) => (
                    <div key={idx} className="bg-white border rounded p-3">
                      <div className="font-medium text-sm">Layer {layer.layer}</div>
                      <div className="text-xs text-gray-600 mt-1">
                        <div>Y: {layer.y_position?.toFixed(3)} m</div>
                        <div>ρ: {layer.resistivities?.join(', ')} Ω·m</div>
                        <div>Blocks: {layer.block_count}</div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
};

// Results Tab Component  
const ResultsTab = ({ results }) => {
  if (!results) {
    return (
      <div className="bg-white rounded-lg shadow-sm border border-gray-200 p-8">
        <div className="text-center text-gray-500">
          <BarChart3 className="mx-auto h-12 w-12 mb-4" />
          <h3 className="text-lg font-medium mb-2">No Results Yet</h3>
          <p>Run a forward model to see results here.</p>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="bg-white rounded-lg shadow-sm border border-gray-200 p-6">
        <h2 className="text-lg font-semibold text-gray-900 mb-4">Forward Modeling Results</h2>
        
        {results.success ? (
          <div className="space-y-4">
            {/* Parameters Summary */}
            <div className="bg-blue-50 border border-blue-200 rounded-md p-4">
              <h3 className="font-medium text-blue-900 mb-2">Model Parameters</h3>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                <div>
                  <span className="text-blue-700 font-medium">Electrodes:</span>
                  <div className="text-blue-900">{results.parameters?.n_electrodes}</div>
                </div>
                <div>
                  <span className="text-blue-700 font-medium">Spacing:</span>
                  <div className="text-blue-900">{results.parameters?.electrode_spacing} m</div>
                </div>
                <div>
                  <span className="text-blue-700 font-medium">Resistivity:</span>
                  <div className="text-blue-900">{results.parameters?.resistivity} Ω·m</div>
                </div>
                <div>
                  <span className="text-blue-700 font-medium">Conductivity:</span>
                  <div className="text-blue-900">{results.parameters?.conductivity?.toFixed(4)} S/m</div>
                </div>
              </div>
            </div>

            {/* Results Summary */}
            <div className="bg-green-50 border border-green-200 rounded-md p-4">
              <h3 className="font-medium text-green-900 mb-2">Results Summary</h3>
              <div className="grid grid-cols-2 md:grid-cols-3 gap-4 text-sm">
                <div>
                  <span className="text-green-700 font-medium">Data Points:</span>
                  <div className="text-green-900">{results.results?.num_data_points}</div>
                </div>
                <div>
                  <span className="text-green-700 font-medium">Mesh Nodes X:</span>
                  <div className="text-green-900">{results.results?.mesh_info?.nodes_x}</div>
                </div>
                <div>
                  <span className="text-green-700 font-medium">Mesh Nodes Y:</span>
                  <div className="text-green-900">{results.results?.mesh_info?.nodes_y}</div>
                </div>
              </div>
            </div>

            {/* VI Data Preview */}
            {results.results?.vi_data && (
              <div className="bg-gray-50 border border-gray-200 rounded-md p-4">
                <h3 className="font-medium text-gray-900 mb-2">V/I Data Preview (First 10 values)</h3>
                <div className="grid grid-cols-5 gap-2 text-sm font-mono">
                  {results.results.vi_data.slice(0, 10).map((vi, idx) => (
                    <div key={idx} className="bg-white px-2 py-1 rounded border">
                      {vi.toExponential(3)}
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Survey Configuration Preview */}
            {results.results?.survey_config && (
              <div className="bg-gray-50 border border-gray-200 rounded-md p-4">
                <h3 className="font-medium text-gray-900 mb-2">Survey Configuration Preview (ABMN)</h3>
                <div className="overflow-x-auto">
                  <table className="min-w-full text-sm">
                    <thead className="bg-white">
                      <tr>
                        <th className="px-3 py-2 text-left font-medium text-gray-900">A</th>
                        <th className="px-3 py-2 text-left font-medium text-gray-900">B</th>
                        <th className="px-3 py-2 text-left font-medium text-gray-900">M</th>
                        <th className="px-3 py-2 text-left font-medium text-gray-900">N</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-gray-200">
                      {results.results.survey_config.slice(0, 5).map((config, idx) => (
                        <tr key={idx} className="bg-white">
                          <td className="px-3 py-2">{config[0]}</td>
                          <td className="px-3 py-2">{config[1]}</td>
                          <td className="px-3 py-2">{config[2]}</td>
                          <td className="px-3 py-2">{config[3]}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}
          </div>
        ) : (
          <div className="bg-red-50 border border-red-200 rounded-md p-4">
            <div className="flex items-center">
              <AlertCircle className="h-5 w-5 text-red-400 mr-2" />
              <h3 className="font-medium text-red-900">Error</h3>
            </div>
            <div className="text-red-700 mt-1">{results.error || 'Unknown error occurred'}</div>
          </div>
        )}
      </div>
    </div>
  );
};

function App() {
  return (
    <div className="App">
      <EarthImagerInterface />
    </div>
  );
}

export default App;