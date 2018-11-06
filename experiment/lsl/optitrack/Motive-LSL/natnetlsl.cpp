#include "natnetlsl.h"
#include <lsl_cpp.h>
#include <winsock.h>
#include <conio.h>
#include <NatNetCAPI.h>

#define dprintf(...) {} /**/

using namespace lsl;

extern sServerDescription g_serverDescription;
extern int g_analogSamplesPerMocapFrame;
extern sNatNetClientConnectParams g_connectParams;
extern NatNetClient* g_pClient;

static stream_outlet* outlet_m = 0;
static stream_outlet* outlet_r = 0;

void natnetlsl_init()
{
	if (outlet_m || outlet_m) {
		printf("resetting LSL outlets\n");
		delete outlet_m;
		delete outlet_r;
	}
	printf("setting up LSL outlets:\n");
	int szmarkers = 15;
	int szrigids = 16;

	char sProgramVersion[255];
	sprintf(sProgramVersion,"NatNet-LSL (build %s %s)", __DATE__, __TIME__);

	char sManufacturer[] = "RU Biophysics";
	char sApplicationVersion[255];
	char sNatNetVersion[255];
	char sClientIP[255];
	char sServerIP[255];
	char sServerName[255];
	char sMocapFrameRate[255];
	char sAnalogSamplesPerMocapFrame[255];

	sprintf(sApplicationVersion, "%s (ver. %d.%d.%d.%d)", g_serverDescription.szHostApp, g_serverDescription.HostAppVersion[0],
		g_serverDescription.HostAppVersion[1], g_serverDescription.HostAppVersion[2], g_serverDescription.HostAppVersion[3]);
	sprintf(sNatNetVersion, "%d.%d.%d.%d", g_serverDescription.NatNetVersion[0], g_serverDescription.NatNetVersion[1],
		g_serverDescription.NatNetVersion[2], g_serverDescription.NatNetVersion[3]);
	sprintf(sClientIP, "%s", g_connectParams.localAddress);
	sprintf(sServerIP, "%s", g_connectParams.serverAddress);
	sprintf(sServerName, "%s", g_serverDescription.szHostComputerName);

	// get mocap frame rate
	void* pResult;
	int nBytes = 0;
	int ret = g_pClient->SendMessageAndWait("FrameRate", &pResult, &nBytes);
	double fRate = lsl::IRREGULAR_RATE;
	if (ret == ErrorCode_OK)
	{
		fRate = *((float*)pResult);
		sprintf(sMocapFrameRate, "%3.2f", fRate);
	}
	else
		sprintf(sMocapFrameRate,"0");

	// get # of analog samples per mocap frame of data
	ret = g_pClient->SendMessageAndWait("AnalogSamplesPerMocapFrame", &pResult, &nBytes);
	if (ret == ErrorCode_OK)
	{
		g_analogSamplesPerMocapFrame = *((int*)pResult);
		sprintf(sAnalogSamplesPerMocapFrame,"%d", g_analogSamplesPerMocapFrame);
	}
	else
		sprintf(sAnalogSamplesPerMocapFrame,"0");

	char hostname[255];
	gethostname(hostname, 255);

	char info_type[255];
	sprintf(info_type, "OptiTrack Mocap @ %s", hostname);

	
	char info_name_labeled_markers[255];
	sprintf(info_name_labeled_markers, "Labeled Markers");

	printf("name=\"%s\"\ttype=\"%s\"\n", info_name_labeled_markers, info_type);
	stream_info info_m(info_name_labeled_markers, info_type, szmarkers, fRate, lsl::channel_format_t(lsl::cft_double64));

	info_m.desc()
		.append_child_value("Version", sProgramVersion)
		.append_child_value("Manufacturer", sManufacturer)
		.append_child_value("MocapVersion", sApplicationVersion)
		.append_child_value("NatNetVersion", sNatNetVersion)
		.append_child_value("ClientIP", sClientIP)
		.append_child_value("ServerIP", sServerIP)
		.append_child_value("ServerName", sServerName);
		//GW: We do not supply the FrameRate here, because it may be changed during an experiment,
		// and we do not want to be forced to restart all lsl streams and have incorrect information if
		// we forget...
		//.append_child_value("MocapFrameRate",sMocapFrameRate)
		//.append_child_value("AnalogSamplesPerMocapFrame",sAnalogSamplesPerMocapFrame);

	const char* labels_m[] = {
		"FrameID", "fTimeStamp", "CameraDataReceivedTimestamp", "TransmitTimestamp", "SoftwareLatency", "TransmitLatency",
		"ModelID", "MarkerID", "Occluded", "PCSolved", "ModelSolved",
		"size", "x", "y", "z"
	};
	lsl::xml_element chns_m = info_m.desc().append_child("fields");
	for (const char* label : labels_m)
	{
		chns_m.append_child("field")
			.append_child_value("label", label);
	}
	outlet_m = new stream_outlet(info_m);

	char info_name_rigid_bodies[255];
	sprintf(info_name_rigid_bodies, "Rigid Bodies");
	printf("name=\"%s\"\ttype=\"%s\"\n", info_name_rigid_bodies, info_type);

	stream_info info_r(info_name_rigid_bodies, info_type, szrigids, fRate, lsl::channel_format_t(lsl::cft_double64));
	info_r.desc()
		.append_child_value("Version", sProgramVersion)
		.append_child_value("Manufacturer", sManufacturer)
		.append_child_value("MocapVersion", sApplicationVersion)
		.append_child_value("NatNetVersion", sNatNetVersion)
		.append_child_value("ClientIP", sClientIP)
		.append_child_value("ServerIP", sServerIP)
		.append_child_value("ServerName", sServerName);
		//.append_child_value("MocapFrameRate", sMocapFrameRate)
		//.append_child_value("AnalogSamplesPerMocapFrame", sAnalogSamplesPerMocapFrame);

	const char* labels_r[] = {
		"FrameID", "fTimeStamp", "CameraDataReceivedTimestamp", "TransmitTimestamp", "SoftwareLatency", "TransmitLatency",
		"BodyID", "Error", "Valid",
		"x", "y", "z", "qx", "qy", "qz", "qw"
	};
	lsl::xml_element chns_r = info_r.desc().append_child("fields");
	for (const char* label : labels_r)
	{
		chns_r.append_child("field")
			.append_child_value("label", label);
	}
	outlet_r = new stream_outlet(info_r);
}

void natnetlsl_write(sFrameOfMocapData* data, NatNetClient* pClient)
{
	dprintf("\nnatnetlsl_write()\n\n");
	// Software latency here is defined as the span of time between:
	//   a) The reception of a complete group of 2D frames from the camera system (CameraDataReceivedTimestamp)
	// and
	//   b) The time immediately prior to the NatNet frame being transmitted over the network (TransmitTimestamp)
	//
	// This figure may appear slightly higher than the "software latency" reported in the Motive user interface,
	// because it additionally includes the time spent preparing to stream the data via NatNet.
	const uint64_t softwareLatencyHostTicks = data->TransmitTimestamp - data->CameraDataReceivedTimestamp;
	const double softwareLatency = softwareLatencyHostTicks / static_cast<double>(g_serverDescription.HighResClockFrequency);
	const double CameraDataReceivedTimestamp = (data->CameraDataReceivedTimestamp) / static_cast<double>(g_serverDescription.HighResClockFrequency);
	const double TransmitTimestamp = (data->TransmitTimestamp) / static_cast<double>(g_serverDescription.HighResClockFrequency);

	// Transit latency is defined as the span of time between Motive transmitting the frame of data, and its reception by the client (now).
	// The SecondsSinceHostTimestamp method relies on NatNetClient's internal clock synchronization with the server using Cristian's algorithm.
	const double transitLatency = pClient->SecondsSinceHostTimestamp(data->TransmitTimestamp);
	int i = 0;

	dprintf("FrameID : %d\n", data->iFrame);
	dprintf("Timestamp : %3.2lf\n", data->fTimestamp);
	dprintf("Software latency : %.5lf seconds\n", softwareLatency);

	// Only recent versions of the Motive software in combination with ethernet camera systems support system latency measurement.
	// If it's unavailable (for example, with USB camera systems, or during playback), this field will be zero.
	const bool bSystemLatencyAvailable = data->CameraMidExposureTimestamp != 0;

	if (bSystemLatencyAvailable)
	{
		// System latency here is defined as the span of time between:
		//   a) The midpoint of the camera exposure window, and therefore the average age of the photons (CameraMidExposureTimestamp)
		// and
		//   b) The time immediately prior to the NatNet frame being transmitted over the network (TransmitTimestamp)
		const uint64_t systemLatencyHostTicks = data->TransmitTimestamp - data->CameraMidExposureTimestamp;
		const double systemLatencyMillisec = (systemLatencyHostTicks * 1000) / static_cast<double>(g_serverDescription.HighResClockFrequency);

		// Client latency is defined as the sum of system latency and the transit time taken to relay the data to the NatNet client.
		// This is the all-inclusive measurement (photons to client processing).
		const double clientLatencyMillisec = pClient->SecondsSinceHostTimestamp(data->CameraMidExposureTimestamp) * 1000.0;

		// You could equivalently do the following (not accounting for time elapsed since we calculated transit latency above):
		//const double clientLatencyMillisec = systemLatencyMillisec + transitLatencyMillisec;

		dprintf("System latency : %.2lf milliseconds\n", systemLatencyMillisec);
		//dprintf("Total client latency : %.2lf milliseconds (transit time +%.2lf ms)\n", clientLatencyMillisec, transitLatencyMillisec);
	}
	else
	{
		dprintf("Transit latency : %.2lf milliseconds\n", transitLatency);
	}

	// FrameOfMocapData params
	bool bIsRecording = ((data->params & 0x01) != 0);
	bool bTrackedModelsChanged = ((data->params & 0x02) != 0);
	if (bIsRecording)
		dprintf("RECORDING\n");
	if (bTrackedModelsChanged)
		dprintf("Models Changed.\n");


	// timecode - for systems with an eSync and SMPTE timecode generator - decode to values
	int hour, minute, second, frame, subframe;
	NatNet_DecodeTimecode(data->Timecode, data->TimecodeSubframe, &hour, &minute, &second, &frame, &subframe);
	// decode to friendly string
	char szTimecode[128] = "";
	NatNet_TimecodeStringify(data->Timecode, data->TimecodeSubframe, szTimecode, 128);
	dprintf("Timecode : %s\n", szTimecode);

	// Rigid Bodies
	dprintf("Rigid Bodies [Count=%d]\n", data->nRigidBodies);
	for (i = 0; i < data->nRigidBodies; i++)
	{
		// params
		// 0x01 : bool, rigid body was successfully tracked in this frame
		bool bTrackingValid = data->RigidBodies[i].params & 0x01;
		// IMPORTANT: lsldata, szrigids and the xml information must be matching !!!!!
		double lsldata[] = { double(data->iFrame), data->fTimestamp, CameraDataReceivedTimestamp, TransmitTimestamp, softwareLatency, transitLatency,
			double(data->RigidBodies[i].ID), data->RigidBodies[i].MeanError, double(bTrackingValid),
			data->RigidBodies[i].x,	data->RigidBodies[i].y,	data->RigidBodies[i].z,
			data->RigidBodies[i].qx,data->RigidBodies[i].qy,data->RigidBodies[i].qz,data->RigidBodies[i].qw
		};
		if (outlet_r)
			outlet_r->push_sample(lsldata);

		dprintf("Rigid Body [ID=%d  Error=%3.2f  Valid=%d]\n", data->RigidBodies[i].ID, data->RigidBodies[i].MeanError, bTrackingValid);
		dprintf("\tx\ty\tz\tqx\tqy\tqz\tqw\n");
		dprintf("\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n",
			data->RigidBodies[i].x,
			data->RigidBodies[i].y,
			data->RigidBodies[i].z,
			data->RigidBodies[i].qx,
			data->RigidBodies[i].qy,
			data->RigidBodies[i].qz,
			data->RigidBodies[i].qw);
	}

	// Skeletons
	dprintf("Skeletons [Count=%d]\n", data->nSkeletons);
	for (i = 0; i < data->nSkeletons; i++)
	{
		sSkeletonData skData = data->Skeletons[i];
		dprintf("Skeleton [ID=%d  Bone count=%d]\n", skData.skeletonID, skData.nRigidBodies);
		for (int j = 0; j< skData.nRigidBodies; j++)
		{
			sRigidBodyData rbData = skData.RigidBodyData[j];
			dprintf("Bone %d\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n",
				rbData.ID, rbData.x, rbData.y, rbData.z, rbData.qx, rbData.qy, rbData.qz, rbData.qw);
		}
	}

	// labeled markers - this includes all markers (Active, Passive, and 'unlabeled' (markers with no asset but a PointCloud ID)
	bool bOccluded;     // marker was not visible (occluded) in this frame
	bool bPCSolved;     // reported position provided by point cloud solve
	bool bModelSolved;  // reported position provided by model solve
	bool bHasModel;     // marker has an associated asset in the data stream
	bool bUnlabeled;    // marker is 'unlabeled', but has a point cloud ID that matches Motive PointCloud ID (In Motive 3D View)
	bool bActiveMarker; // marker is an actively labeled LED marker

	dprintf("Markers [Count=%d]\n", data->nLabeledMarkers);
	for (i = 0; i < data->nLabeledMarkers; i++)
	{
		bOccluded = ((data->LabeledMarkers[i].params & 0x01) != 0);
		bPCSolved = ((data->LabeledMarkers[i].params & 0x02) != 0);
		bModelSolved = ((data->LabeledMarkers[i].params & 0x04) != 0);
		bHasModel = ((data->LabeledMarkers[i].params & 0x08) != 0);
		bUnlabeled = ((data->LabeledMarkers[i].params & 0x10) != 0);
		bActiveMarker = ((data->LabeledMarkers[i].params & 0x20) != 0);

		sMarker marker = data->LabeledMarkers[i];

		// Marker ID Scheme:
		// Active Markers:
		//   ID = ActiveID, correlates to RB ActiveLabels list
		// Passive Markers: 
		//   If Asset with Legacy Labels
		//      AssetID 	(Hi Word)
		//      MemberID	(Lo Word)
		//   Else
		//      PointCloud ID
		int modelID, markerID;
		NatNet_DecodeID(marker.ID, &modelID, &markerID);

		char szMarkerType[512];
		if (bActiveMarker)
			strcpy(szMarkerType, "Active");
		else if (bUnlabeled)
			strcpy(szMarkerType, "Unlabeled");
		else
			strcpy(szMarkerType, "Labeled");

		dprintf("%s Marker [ModelID=%d, MarkerID=%d, Occluded=%d, PCSolved=%d, ModelSolved=%d] [size=%3.2f] [pos=%3.2f,%3.2f,%3.2f]\n",
			szMarkerType, modelID, markerID, bOccluded, bPCSolved, bModelSolved, marker.size, marker.x, marker.y, marker.z);
		// IMPORTANT: lsldata, szmarkers and the xml information must be matching !!!!!
		double lsldata[] = { double(data->iFrame), data->fTimestamp, CameraDataReceivedTimestamp, TransmitTimestamp, softwareLatency, transitLatency,
			double(modelID), double(markerID), double(bOccluded), double(bPCSolved), double(bModelSolved),
			marker.size, marker.x, marker.y, marker.z
		};
		if (outlet_m)
			outlet_m->push_sample(lsldata);
	}

	// force plates
	dprintf("Force Plate [Count=%d]\n", data->nForcePlates);
	for (int iPlate = 0; iPlate < data->nForcePlates; iPlate++)
	{
		printf("Force Plate %d\n", data->ForcePlates[iPlate].ID);
		for (int iChannel = 0; iChannel < data->ForcePlates[iPlate].nChannels; iChannel++)
		{
			printf("\tChannel %d:\t", iChannel);
			if (data->ForcePlates[iPlate].ChannelData[iChannel].nFrames == 0)
			{
				printf("\tEmpty Frame\n");
			}
			else if (data->ForcePlates[iPlate].ChannelData[iChannel].nFrames != g_analogSamplesPerMocapFrame)
			{
				printf("\tPartial Frame [Expected:%d   Actual:%d]\n", g_analogSamplesPerMocapFrame, data->ForcePlates[iPlate].ChannelData[iChannel].nFrames);
			}
			for (int iSample = 0; iSample < data->ForcePlates[iPlate].ChannelData[iChannel].nFrames; iSample++)
				printf("%3.2f\t", data->ForcePlates[iPlate].ChannelData[iChannel].Values[iSample]);
			printf("\n");
		}
	}

	// devices
	dprintf("Device [Count=%d]\n", data->nDevices);
	for (int iDevice = 0; iDevice < data->nDevices; iDevice++)
	{
		printf("Device %d\n", data->Devices[iDevice].ID);
		for (int iChannel = 0; iChannel < data->Devices[iDevice].nChannels; iChannel++)
		{
			printf("\tChannel %d:\t", iChannel);
			if (data->Devices[iDevice].ChannelData[iChannel].nFrames == 0)
			{
				printf("\tEmpty Frame\n");
			}
			else if (data->Devices[iDevice].ChannelData[iChannel].nFrames != g_analogSamplesPerMocapFrame)
			{
				printf("\tPartial Frame [Expected:%d   Actual:%d]\n", g_analogSamplesPerMocapFrame, data->Devices[iDevice].ChannelData[iChannel].nFrames);
			}
			for (int iSample = 0; iSample < data->Devices[iDevice].ChannelData[iChannel].nFrames; iSample++)
				printf("%3.2f\t", data->Devices[iDevice].ChannelData[iChannel].Values[iSample]);
			printf("\n");
		}
	}
}

void natnetlsl_done()
{
	delete outlet_m;
	delete outlet_r;
}


